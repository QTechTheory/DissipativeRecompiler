
#include "trotterevolver.h"
#include "hamiltonianloader.h"
#include "circuitloader.h"

#include <QuEST.h>


#define TROT_MODE_FIXED 0
#define TROT_MODE_REVERSE 1
#define TROT_MODE_SHUFFLE 2


void getTrotterCycleParams(double* params, Hamiltonian hamil, double duration, int numCycles) {
    
    // e^(-i duration coeff X) = Rx( 2 * duration * coeff / numCycles )^numCycles
    for (int t=0; t < hamil.numTerms; t++)
        params[t] = 2 * duration * hamil.termCoeffs[t] / (double) numCycles;
}

/**
 * Requires that the terms of hamil correspond to the gates of trotterCycle
 */
void produceTrotterState(
    QubitRegister initState, QubitRegister outState, 
    Hamiltonian hamil, Circuit trotterCycle,
    double* cycleParams, double duration, int numCycles, int mode
) {
    if (trotterCycle.numGates != hamil.numTerms) {
        printf("ERROR! Given trotter cycle circuit has a different number of gates than Hamiltonian terms!\n");
        exit(1);
    }
    
    if (mode != TROT_MODE_FIXED && mode != TROT_MODE_REVERSE && mode != TROT_MODE_SHUFFLE) {
        printf(
            "ERROR! Trotter mode %d unrecognised! Must be 0 (fixed), "
            "1 (every 2nd cycle reversed) or "
            "2 (each cycle is randomly shuffled)! Exiting...\n", mode);
        exit(1);
    }
    
    cloneQubitRegister(outState, initState);
    getTrotterCycleParams(cycleParams, hamil, duration, numCycles);
    
    for (int c=0; c < numCycles; c++) {
        
        // apply set order of terms
        if (mode == TROT_MODE_FIXED)
            applyCircuit(trotterCycle, cycleParams, outState);
        
        // shuffle each Trotter cycle
        if (mode == TROT_MODE_SHUFFLE)
            applyRandomisedCircuit(trotterCycle, cycleParams, outState);
        
        // reverse every second cycle
        if (mode == TROT_MODE_REVERSE) {
            if (c % 2 == 0)
                applyCircuit(trotterCycle, cycleParams, outState);
            else
                applyReversedCircuit(trotterCycle, cycleParams, outState);
        }
    }
}