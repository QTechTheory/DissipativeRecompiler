
#include "paramevolver.h"
#include "circuitloader.h"

#include <QuEST.h>
#include <stdio.h>
#include <stdlib.h>


#define REAL_SOLVER_METHOD 1            // 0 for TSVD, 1 for Tikhonov
#define IMAG_SOLVER_METHOD 0            // 0 for TSVD, 1 for Tikhonov
#define TSVD_TOLERANCE 1E-5
#define TIKHONOV_PARAM_SEARCH_SIZE 3    // must be >= 3

#define DERIV_ORDER 4
#define DERIV_STEP_SIZE 1E-5

double FIRST_DERIV_FINITE_DIFFERENCE_COEFFS[4][4] = {
	{1/2.0},
	{2/3.0, -1/12.0},
	{3/4.0, -3/20.0, 1/60.0},
	{4/5.0, -1/5.0, 4/105.0, -1/280.0}
};


ParamEvolEnv initParamEvolEnv(int numQubits, int numParams, QuESTEnv qEnv) {
    
    ParamEvolEnv evEnv;
    evEnv.numParams = numParams;
    evEnv.hamilWavef = createQubitRegister(numQubits, qEnv);
    evEnv.copyWavef = createQubitRegister(numQubits, qEnv);
    
    evEnv.derivWavefs = malloc(numParams * sizeof *evEnv.derivWavefs);
    for (int p=0; p < numParams; p++)
        evEnv.derivWavefs[p] = createQubitRegister(numQubits, qEnv);
    
    evEnv.varMatrix = gsl_matrix_alloc(numParams, numParams);
    evEnv.varVector = gsl_vector_alloc(numParams);
    evEnv.paramChange = gsl_vector_alloc(numParams);
    evEnv.linearSolverSpace = gsl_multifit_linear_alloc(numParams, numParams);
    
    // TSVD
    evEnv.tsvdCovarMatrix = gsl_matrix_alloc(numParams, numParams);
    
    // Tikhonov
	evEnv.tikhonovParamSamples = gsl_vector_alloc(TIKHONOV_PARAM_SEARCH_SIZE);
	evEnv.tikhonovParamRho = gsl_vector_alloc(TIKHONOV_PARAM_SEARCH_SIZE);
	evEnv.tikhonovParamEta = gsl_vector_alloc(TIKHONOV_PARAM_SEARCH_SIZE);
	evEnv.tikhonovVecL = gsl_vector_alloc(numParams);
	for (int i=0; i < numParams; i++)
		gsl_vector_set(evEnv.tikhonovVecL, i, 1);

    return evEnv;
}

void closeParamEvolEnv(ParamEvolEnv evEnv, QuESTEnv qEnv) {
    
    destroyQubitRegister(evEnv.hamilWavef, qEnv);
    destroyQubitRegister(evEnv.copyWavef, qEnv);
    
    for (int p=0; p < evEnv.numParams; p++)
        destroyQubitRegister(evEnv.derivWavefs[p], qEnv);
    free(evEnv.derivWavefs);
    
    gsl_matrix_free(evEnv.varMatrix);
    gsl_vector_free(evEnv.varVector);
    gsl_vector_free(evEnv.paramChange);
    gsl_multifit_linear_free(evEnv.linearSolverSpace);
    
    // TSVD
    gsl_matrix_free(evEnv.tsvdCovarMatrix);
    
    // Tikhonov
	gsl_vector_free(evEnv.tikhonovParamSamples);
	gsl_vector_free(evEnv.tikhonovParamRho);
	gsl_vector_free(evEnv.tikhonovParamEta);
	gsl_vector_free(evEnv.tikhonovVecL);
}

void populateDerivs(
    ParamEvolEnv evEnv, QubitRegister initState, Circuit ansatz, double* params
) {
    
    for (int p=0; p < evEnv.numParams; p++) {
        
        // clear deriv
        QubitRegister deriv = evEnv.derivWavefs[p];
        for (long long int i=0; i < deriv.numAmpsTotal; i++) {
            deriv.stateVec.real[i] = 0;
            deriv.stateVec.imag[i] = 0;
        }
        
    	// approx deriv with finite difference
    	double* coeffs = FIRST_DERIV_FINITE_DIFFERENCE_COEFFS[DERIV_ORDER - 1];
        double origParam = params[p];
        
    	// repeatly add c*psi(p+ndp) - c*psi(p-ndp) to deriv
    	for (int step=1; step <= DERIV_ORDER; step++) {
            for (int sign = -1; sign <= 1; sign+=2) {
                params[p] = origParam + sign*step*DERIV_STEP_SIZE;
                
                cloneQubitRegister(evEnv.copyWavef, initState);
                applyCircuit(ansatz, params, evEnv.copyWavef);
                
                for (long long int i=0; i < deriv.numAmpsTotal; i++) {
                    deriv.stateVec.real[i] += (sign * coeffs[step-1] * 
                        getRealAmpEl(evEnv.copyWavef, i));
                    deriv.stateVec.imag[i] += (sign * coeffs[step-1] * 
                        getImagAmpEl(evEnv.copyWavef, i));
                }
            }
        }
        
        // divide by the step size
        for (long long int i=0; i < deriv.numAmpsTotal; i++) {
            deriv.stateVec.real[i] /= DERIV_STEP_SIZE;
            deriv.stateVec.imag[i] /= DERIV_STEP_SIZE;
        }
        
        // restore the parameter
        params[p] = origParam;   
    }
}

Complex getConj(Complex a) {
    Complex conj = a;
    conj.imag *= -1;
    return conj;
}

void populateMatrices(
    ParamEvolEnv evEnv, QubitRegister initState, Circuit ansatz, double* params, 
    Hamiltonian hamil, int forRealTime
) {
    // compute |dpsi/dp>
    populateDerivs(evEnv, initState, ansatz, params);
    
    // compute H|psi>
    cloneQubitRegister(evEnv.copyWavef, initState);
    applyCircuit(ansatz, params, evEnv.copyWavef);
    applyHamiltonian(hamil, evEnv.copyWavef, evEnv.hamilWavef);
    
    // populate <dpsi/dp_i | dpsi/dp_j>
    for (int i=0; i < evEnv.numParams; i++) {
        for (int j=0; j < evEnv.numParams; j++) {
            
            // real-time takes imag component, which is zero on diagonals
            if (i == j && forRealTime) {
                gsl_matrix_set(evEnv.varMatrix, i, j, 0);
                // (imag-time norm depends on whether gate is controlled, so no shortcuts!)
            }
            // conjugate of previously calculated inner product (matr is conj-symmetric)
            else if (i > j) {
                double prev = gsl_matrix_get(evEnv.varMatrix, j, i);
                double comp = (forRealTime)? -prev : prev;
                gsl_matrix_set(evEnv.varMatrix, i, j, comp);
            }
            // get component of deriv inner product
            else {
                Complex prod = calcInnerProduct(evEnv.derivWavefs[i], evEnv.derivWavefs[j]);
                double comp = (forRealTime)? -prod.imag : prod.real;
                gsl_matrix_set(evEnv.varMatrix, i, j, comp);
            }
        }
    }
    
    // populate <dpsi/dp_i|H|psi>
    for (int i=0; i < evEnv.numParams; i++) {
        Complex prod = calcInnerProduct(evEnv.derivWavefs[i], evEnv.hamilWavef);
        double comp = (forRealTime)? prod.real : -prod.real;
        gsl_vector_set(evEnv.varVector, i, comp);
    }
}

double solveViaTSVD(ParamEvolEnv evEnv) {
    
    double residSum; 
    size_t singValsKept;
    gsl_multifit_linear_tsvd(
        evEnv.varMatrix, evEnv.varVector, TSVD_TOLERANCE, 
        evEnv.paramChange, evEnv.tsvdCovarMatrix, &residSum, &singValsKept, 
        evEnv.linearSolverSpace);
    
    return residSum;
}

double solveViaTikhonov(ParamEvolEnv evEnv) {

	// compute the SVD in the workspace (needed for L-curve and solve)
    gsl_multifit_linear_svd(evEnv.varMatrix, evEnv.linearSolverSpace);
    
	// sample the system under different regularisation params (build L-curve)
    gsl_multifit_linear_lcurve(
    	evEnv.varVector, evEnv.tikhonovParamSamples,
    	evEnv.tikhonovParamRho, evEnv.tikhonovParamEta, evEnv.linearSolverSpace);
            
	// choose the best regularisation param (hardcode 0.02 performs ok)
    size_t tikhonovParamIndex;
    gsl_multifit_linear_lcorner(
    	evEnv.tikhonovParamRho, evEnv.tikhonovParamEta, &tikhonovParamIndex);
    double tikhonovParam = gsl_vector_get(evEnv.tikhonovParamSamples, tikhonovParamIndex);
    
	// the error ||varVector - varMatrix paramChange|| can be monitored
	double residualNorm, paramChangeNorm;
    
	// perform Tikhonov regularisation (where L = identity)
	gsl_multifit_linear_solve(
		tikhonovParam, evEnv.varMatrix, evEnv.varVector, evEnv.paramChange,
		&residualNorm, &paramChangeNorm, evEnv.linearSolverSpace);
        
    return residualNorm;
}

double evolveParamsInner(
    ParamEvolEnv evEnv, QubitRegister initState, Circuit ansatz, double* params,
    Hamiltonian hamil, int inRealTime, double timeStep, 
    double (*solveMatrices)(ParamEvolEnv)
) {
    
    // safety first!
    if (evEnv.numParams != ansatz.numGates) {
        printf(
            "ERROR! ParamEvolEnv has a different number of params than the ansatz has gates! "
            "Exiting...\n");
        exit(1);
    }
    
    // populate varMatrix and varVector
    populateMatrices(evEnv, initState, ansatz, params, hamil, inRealTime);
    
    // populate paramChange
    double error = solveMatrices(evEnv);

    // update the params
    for (int p=0; p < evEnv.numParams; p++)
        params[p] += timeStep * gsl_vector_get(evEnv.paramChange, p);
    
    return error;
}

double evolveParamsReal(
    double* params, QubitRegister initState, Circuit ansatz, 
    Hamiltonian hamil, double timeStep, ParamEvolEnv evEnv
) {
    return evolveParamsInner(
        evEnv, initState, ansatz, params, hamil, 1, timeStep, 
        (REAL_SOLVER_METHOD == 0)? solveViaTSVD : solveViaTikhonov);
}
double evolveParamsImag(
    double* params, QubitRegister initState, Circuit ansatz, 
    Hamiltonian hamil, double timeStep, ParamEvolEnv evEnv
) {
    return evolveParamsInner(
        evEnv, initState, ansatz, params, hamil, 0, timeStep,
        (IMAG_SOLVER_METHOD == 0)? solveViaTSVD : solveViaTikhonov);
}
