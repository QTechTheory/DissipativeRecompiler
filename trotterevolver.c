
#include "trotterevolver.h"
#include "hamiltonianloader.h"
#include "circuitloader.h"

#include <QuEST.h>


#define TROT_MODE_FIXED 0
#define TROT_MODE_REVERSE 1
#define TROT_MODE_SHUFFLE 2
#define TROT_MODE_SHUFFLE_REVERSE 3


void getTrotterCycleParams(double* params, Hamiltonian hamil, double duration, int numCycles) {
    
    // e^(-i duration coeff X) = Rx( 2 * duration * coeff / numCycles )^numCycles
    for (int t=0; t < hamil.numTerms; t++)
        params[t] = 2 * duration * hamil.termCoeffs[t] / (double) numCycles;
}


void shuffle(int* indices, int len) {
        
    int indIncluded[len];
    for (int i=0; i < len; i++)
        indIncluded[i] = 0;
    
    int indicesCopy[len];
    for (int i=0; i < len; i++) {
        
        int j = -1;
        while (j == -1 || indIncluded[j])
            j = (rand() * (long long int) len) / RAND_MAX;

        indicesCopy[i] = indices[j];
        indIncluded[j] = 1;
    }
    
    for (int i=0; i < len; i++)
        indices[i] = indicesCopy[i];
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
    
    if (mode != TROT_MODE_FIXED && 
        mode != TROT_MODE_REVERSE && 
        mode != TROT_MODE_SHUFFLE && 
        mode != TROT_MODE_SHUFFLE_REVERSE) {
        printf(
            "ERROR! Trotter mode %d unrecognised! Must be 0 (fixed), "
            "1 (every 2nd cycle reversed) or "
            "2 (each cycle is randomly shuffled)! Exiting...\n", mode);
        exit(1);
    }
    
    // used only by TROT_MODE_SHUFFLE_REVERSE
    int gateOrdering[trotterCycle.numGates];
    for (int i=0; i < trotterCycle.numGates; i++)
        gateOrdering[i] = i;
    
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
        
        // shuffle the circuit after alternatingly reversing
        if (mode == TROT_MODE_SHUFFLE_REVERSE) {
            if (c % 2 == 0)
                applyOrderedCircuit(trotterCycle, cycleParams, gateOrdering, outState);
            else {
                applyReversedOrderedCircuit(trotterCycle, cycleParams, gateOrdering, outState);
                shuffle(gateOrdering, trotterCycle.numGates);
            }
        }
    }
}