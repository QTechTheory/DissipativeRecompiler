#ifndef TROTTER_EVOLVER_H_
#define TROTTER_EVOLVER_H_

#include <QuEST.h>

#include "hamiltonianloader.h"
#include "circuitloader.h"

void produceTrotterState(
    QubitRegister initState, QubitRegister outState, 
    Hamiltonian hamil, Circuit trotterCycle,
    double* cycleParams, double duration, int numCycles,
    int mode
);

#endif // TRUE_EVOLVER_H_