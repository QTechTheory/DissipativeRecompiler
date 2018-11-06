#ifndef PARAM_EVOLVER_H_
#define PARAM_EVOLVER_H_

#include <QuEST.h>
#include <math.h>
#include <complex.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>

#include "circuitloader.h"
#include "hamiltonianloader.h"

typedef struct {
    
    int numParams;
    QubitRegister* derivWavefs;
    QubitRegister hamilWavef;
    QubitRegister copyWavef;
    
    gsl_matrix* varMatrix;
    gsl_vector* varVector;
    gsl_vector* paramChange;
    
    // for TVSD
    gsl_multifit_linear_workspace *linearSolverSpace;
    gsl_matrix *tsvdCovarMatrix;
    
    // for Tikhonov
	gsl_vector *tikhonovParamSamples;
	gsl_vector *tikhonovParamRho;
	gsl_vector *tikhonovParamEta;
	gsl_vector *tikhonovVecL;
    
} ParamEvolEnv;

ParamEvolEnv initParamEvolEnv(int numQubits, int numParams, QuESTEnv qEnv);

double evolveParamsReal(
    double* params, QubitRegister initState, Circuit ansatz, Hamiltonian hamil, 
    double timeStep, ParamEvolEnv evEnv
);
    
double evolveParamsImag(
    double* params, QubitRegister initState, Circuit ansatz, Hamiltonian hamil, 
    double timeStep, ParamEvolEnv evEnv
);
    
void closeParamEvolEnv(ParamEvolEnv evEnv, QuESTEnv qEnv);



/* exposed only for gate eliminator */

void populateMatrices(
    ParamEvolEnv evEnv, QubitRegister initState, Circuit ansatz, double* params, 
    Hamiltonian hamil, int forRealTime
);
    
double solveViaTikhonov(ParamEvolEnv evEnv);

double solveViaTSVD(ParamEvolEnv evEnv);



#endif // PARAM_EVOLVER_H_