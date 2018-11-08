#include <QuEST.h>
#include <stdio.h>
#include <stdlib.h>

#include "circuitloader.h"
#include "hamiltonianloader.h"
#include "paramevolver.h"
#include "trotterevolver.h"
#include "trueevolver.h"
#include "utilities.h"
#include "mmaformatter.h"

#define REAL_TIME_STEP .0025
#define IMAG_TIME_STEP .01
#define REAL_NUM_ITERS 700
#define IMAG_NUM_ITERS 500

int main(int narg, char *varg[]) {
    


    /* load the system info */
    
    QuESTEnv questEnv = initQuESTEnv();
    Circuit oldAnsatz = loadCircuit("data/ansatz0.txt");
    Circuit newAnsatz = loadCircuit("data/ansatz1.txt");
    double* oldParams = loadParams("data/params0_init.txt", oldAnsatz.numGates);
    double* newParams = loadParams("data/params1_init.txt", newAnsatz.numGates);
    Hamiltonian simHamil = loadHamiltonian("data/hamil_sim.txt", questEnv);
    Hamiltonian recHamil = loadHamiltonian("data/hamil_rec.txt", questEnv);
    ParamEvolEnv oldParamEvolver = initParamEvolEnv(oldAnsatz.numQubits, oldAnsatz.numGates, questEnv);
    ParamEvolEnv newParamEvolver = initParamEvolEnv(newAnsatz.numQubits, newAnsatz.numGates, questEnv);
    TrueEvolEnv trueEvolver = initTrueEvolEnv(simHamil, REAL_TIME_STEP, questEnv);

    // display system info
    printf("\nOld ansatz:\n\n");
    printCircuit(oldAnsatz);
    printf("\nNew ansatz:\n\n");
    printCircuit(newAnsatz);
    printf("\nSimulation Hamiltonian:\n\n");
    printHamiltonian(simHamil);
    printf("\nRecompilation Hamiltonian:\n\n");
    printHamiltonian(recHamil);
    
    

    /* allocate state-vectors */
    
    // |in> = |++++++1>
    QubitRegister initState = createQubitRegister(oldAnsatz.numQubits, questEnv);
    initStatePlus(initState);
    hadamard(initState, 0);
    sigmaX(initState,   0);
    
    // |old> = oldAnsatz(oldParams) |in>
    QubitRegister oldState = createQubitRegister(oldAnsatz.numQubits, questEnv);
    cloneQubitRegister(oldState, initState);
    applyCircuit(oldAnsatz, oldParams, oldState);
    
    // |new> = newAnsatz(newParams) |old>
    QubitRegister newState = createQubitRegister(newAnsatz.numQubits, questEnv);
    cloneQubitRegister(newState, oldState);
    applyCircuit(newAnsatz, newParams, newState);
    
    // |true> = timeEvolution(t) |in>
    QubitRegister trueState = createQubitRegister(oldAnsatz.numQubits, questEnv);
    cloneQubitRegister(trueState, initState);
    
    
    
    /* evolve old parameters under real-time Li simulation */
    
    printf("\nSimulating realtime with Li's algorithm and the old ansatz...\n\n");
    
    for (int iter=0; iter < REAL_NUM_ITERS; iter++) {
        
        // update the old parameters
        evolveParamsReal(oldParams, initState, oldAnsatz, simHamil, REAL_TIME_STEP, oldParamEvolver);
        
        // update |old>
        cloneQubitRegister(oldState, initState);
        applyCircuit(oldAnsatz, oldParams, oldState);
        
        // update |true>
        evolveTrueState(trueState, simHamil, trueEvolver);
        
        // monitor |<old|true>|^2
        printf("t=%d, |<old|true>|^2 = %.10f\n", iter, calcFidelity(trueState, oldState));
    }
    

    
    /* evolve new parameters under imaginary-time recompilation */
    
    printf("\nRecompiling the old ansatz into the new ansatz...\n\n");
    
    for (int iter=0; iter < IMAG_NUM_ITERS; iter++) {
        
        // update the new parameters
        evolveParamsImag(newParams, oldState, newAnsatz, recHamil, IMAG_TIME_STEP, newParamEvolver);
        
        // update |new>
        cloneQubitRegister(newState, oldState);
        applyCircuit(newAnsatz, newParams, newState);
        
        // monitor |<old|new>|^2 = |<in|new|old>|^2 and <H>
        printf("t=%d, |<old|new>|^2 = %.5f, <H>-E0 = %.5f\n", iter, 
            calcFidelity(newState, initState), 
            getEnergy(recHamil, newState) - recHamil.eigvals[0]);
    }
    
    
    
    /* free data structures */
    
    free(oldParams);
    free(newParams);
    freeCircuit(oldAnsatz);
    freeCircuit(newAnsatz);
    freeHamiltonian(simHamil, questEnv);
    freeHamiltonian(recHamil, questEnv);
    destroyQubitRegister(oldState, questEnv);
    destroyQubitRegister(newState, questEnv);
    destroyQubitRegister(initState, questEnv);
    destroyQubitRegister(trueState, questEnv);
    closeParamEvolEnv(oldParamEvolver, questEnv);
    closeParamEvolEnv(newParamEvolver, questEnv);
    closeTrueEvolEnv(trueEvolver, questEnv);
    closeQuESTEnv(questEnv); 
    return EXIT_SUCCESS;
}