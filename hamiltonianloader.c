#include "hamiltonianloader.h"
#include "circuitloader.h"

#include <QuEST.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>

/* signatures */

void populateHamiltonianMatrix(Hamiltonian hamil, QuESTEnv env);
double* diagonaliseHamiltonianMatrix(Hamiltonian hamil, QuESTEnv env);

/* definitions */

Hamiltonian loadHamiltonian(char* filename, QuESTEnv env) {
    
    // read in pauli terms
    Hamiltonian hamil;
    hamil.termGates = loadCircuitFromHamilFile(filename,1);
    hamil.numQubits = hamil.termGates.numQubits;
    hamil.numTerms = hamil.termGates.numGates;
    hamil.termCoeffs = malloc(hamil.numTerms * sizeof *(hamil.termCoeffs));
    hamil.copyReg = createQubitRegister(hamil.numQubits, env);
        
    // re-open file
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
		printf("ERROR: file (%s) not found!\n", filename);
        exit(1);
    }
    
    // read in coefficients
    rewind(file);
    for (int t=0; t < hamil.numTerms; t++) {
        fscanf(file, "%lf ", &(hamil.termCoeffs[t]));
        fscanf(file, "%*[^\n]"); // jump to next line
    }
    fclose(file);
    
    // convert pauli terms to matrix
    long long int numAmps = 2LL << (hamil.numQubits - 1);
    hamil.matrix = gsl_matrix_complex_calloc(numAmps, numAmps);
    populateHamiltonianMatrix(hamil, env);
    
    // diagonalise matrix
    hamil.eigvals = diagonaliseHamiltonianMatrix(hamil, env);
    
    return hamil;
}

void printHamiltonian(Hamiltonian hamil) {
    
    printf(
        "The Hamiltonian has %d terms acting on %d qubits. They are...\n",
        hamil.numTerms, hamil.numQubits);
    
    for (int t=0; t < hamil.numTerms; t++) {
        printf("%lf of ", hamil.termCoeffs[t]);
        printGate(hamil.termGates, t);
    }
}

void freeHamiltonian(Hamiltonian hamil, QuESTEnv env) {
    
    free(hamil.eigvals);
    free(hamil.termCoeffs);
    freeCircuit(hamil.termGates);
    destroyQubitRegister(hamil.copyReg, env);
    gsl_matrix_complex_free(hamil.matrix);
}

double getEnergy(Hamiltonian hamil, QubitRegister qureg) {
    
    double energy = 0;
    cloneQubitRegister(hamil.copyReg, qureg);
    
    for (int t=0; t < hamil.numTerms; t++) {
        applyGate(hamil.termGates, t, hamil.copyReg);
    
        double inner = 0;
        for (long long int i=0; i < qureg.numAmpsTotal; i++) {
            double a1 = getRealAmpEl(qureg, i);
            double b1 = - getImagAmpEl(qureg, i);
            double a2 = getRealAmpEl(hamil.copyReg, i);
            double b2 = getImagAmpEl(hamil.copyReg, i);
            inner += a1*a2 - b1*b2;
        }
        
        energy += hamil.termCoeffs[t] * inner;
        cloneQubitRegister(hamil.copyReg, qureg);
    }
    
    return energy;
}

void applyHamiltonian(Hamiltonian hamil, QubitRegister qureg, QubitRegister out) {

    // clear out
    for (long long int i=0; i < out.numAmpsTotal; i++) {
        out.stateVec.real[i] = 0;
        out.stateVec.imag[i] = 0;
    }
    
    // out = hamil qureg
    cloneQubitRegister(hamil.copyReg, qureg);

    //     = sum_t coeff_t term_t qureg
    for (int t=0; t < hamil.numTerms; t++) {
        applyGate(hamil.termGates, t, hamil.copyReg);
        
        double coeff = hamil.termCoeffs[t];
        for (long long int i=0; i < out.numAmpsTotal; i++) {
            
            out.stateVec.real[i] += coeff * hamil.copyReg.stateVec.real[i];
            out.stateVec.imag[i] += coeff * hamil.copyReg.stateVec.imag[i];
        }
        
        cloneQubitRegister(hamil.copyReg, qureg);
    }
}

void populateHamiltonianMatrix(Hamiltonian hamil, QuESTEnv env) {
    
    QubitRegister basisReg = createQubitRegister(hamil.numQubits, env);
    QubitRegister hamilReg = createQubitRegister(hamil.numQubits, env);
    
    // the j-th column of hamil is the result of applying hamil to basis state |j>
    
    long long int numAmps = 2LL << (hamil.numQubits - 1);
        
    for (long long int j=0; j < numAmps; j++) {
        initClassicalState(basisReg, j);
        applyHamiltonian(hamil, basisReg, hamilReg);
        
        for (long long int i=0; i < numAmps; i++) {
            gsl_complex amp = gsl_complex_rect(
                getRealAmpEl(hamilReg,i), getImagAmpEl(hamilReg,i));
            
            gsl_matrix_complex_set(hamil.matrix, i, j, amp);
        }
    }
    
    destroyQubitRegister(basisReg, env);
    destroyQubitRegister(hamilReg, env);
}

double* diagonaliseHamiltonianMatrix(Hamiltonian hamil, QuESTEnv env) {
            
    // prepare GSL data structures
    long long int numAmps = 2LL << (hamil.numQubits - 1);
	gsl_eigen_hermv_workspace* space = gsl_eigen_hermv_alloc(numAmps);
	gsl_vector* eigValsVec = gsl_vector_alloc(numAmps);
	gsl_matrix_complex* eigVecsMatr = gsl_matrix_complex_alloc(numAmps, numAmps); 
    
    // copy hamil matrix (it's damaged by diagonalisation)
    gsl_matrix_complex* hamilCopy = gsl_matrix_complex_alloc(numAmps, numAmps);
    gsl_matrix_complex_memcpy(hamilCopy, hamil.matrix);
        
    // diagonalise H matrix
    int failed = gsl_eigen_hermv(hamil.matrix, eigValsVec, eigVecsMatr, space);
    if (failed) {
        printf("Hamiltonian diagonalisation (through GSL) failed! Exiting...\n");
        exit(1);
    }
    
    // restore damaged hamil matrix
    gsl_matrix_complex_memcpy(hamil.matrix, hamilCopy);
    gsl_matrix_complex_free(hamilCopy);
        
    // sort spectrum by increasing energy
	gsl_eigen_genhermv_sort(eigValsVec, eigVecsMatr, GSL_EIGEN_SORT_VAL_ASC);
    
    // copy from GSL object to pointer
    double* eigvals = malloc(numAmps * sizeof *eigvals);
    for (int i=0; i < numAmps; i++)
        eigvals[i] = gsl_vector_get(eigValsVec, i);    
    
	// free GSL objects
	gsl_eigen_hermv_free(space);
	gsl_vector_free(eigValsVec);
	gsl_matrix_complex_free(eigVecsMatr);
    
    return eigvals;  
}






