#ifndef HAMILTONIAN_LOADER_H_
#define HAMILTONIAN_LOADER_H_

#include <QuEST.h>
#include <gsl/gsl_matrix.h>

#include "circuitloader.h"

typedef struct {
    
    int numQubits;
    int numTerms;
    
    double* termCoeffs;
    Circuit termGates;
    
    gsl_matrix_complex* matrix;
    double* eigvals;
    
    QubitRegister copyReg;
    
} Hamiltonian;

Hamiltonian loadHamiltonian(char* filename, QuESTEnv env);

void printHamiltonian(Hamiltonian hamil);

void applyHamiltonian(Hamiltonian hamil, QubitRegister qureg, QubitRegister out);

void freeHamiltonian(Hamiltonian hamil, QuESTEnv env);

double getEnergy(Hamiltonian hamil, QubitRegister qureg);

#endif // HAMILTONIAN_LOADER_H_