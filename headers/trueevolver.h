#ifndef TRUE_EVOLVER_H_
#define TRUE_EVOLVER_H_

#include <QuEST.h>
#include <gsl/gsl_matrix.h>

#include "hamiltonianloader.h"

typedef struct {
    
    int numQubits;
    long long int dim;
    double timeStep;
    
    gsl_matrix_complex* propogator;
    gsl_vector_complex* copyStateInVec;
    gsl_vector_complex* copyStateOutVec;
    
} TrueEvolEnv;

TrueEvolEnv initTrueEvolEnv(Hamiltonian hamil, double timeStep, QuESTEnv qEnv);

void evolveTrueState(QubitRegister qureg, Hamiltonian hamil, TrueEvolEnv evEnv);

void closeTrueEvolEnv(TrueEvolEnv evEnv, QuESTEnv qEnv);

#endif // TRUE_EVOLVER_H_