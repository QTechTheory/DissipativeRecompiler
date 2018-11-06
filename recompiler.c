#include <QuEST.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

#include "circuitloader.h"
#include "hamiltonianloader.h"
#include "paramevolver.h"
#include "utilities.h"
#include "mmaformatter.h"



#define RECOMP_HAMIL_FN "data/hamil_recomp.txt"

#define IMAG_EVO_TIME_STEP 0.01
#define IMAG_EVO_NUM_ITERS 500
#define IMAG_EVO_INIT_TARGET_PARAM 1E-8

#define OUTPUT_NUM_SIGFIGS 10



int main(int narg, char *varg[]) {
    
    /* collect arguments */

    if (narg != 1 + 9) {
        printf("args:\n"
            "old_ansatz_fn\n"
            "old_params_fn\n"
            "new_ansatz_fn\n"
            "new_init_params_fn\n"
            "true_state_fn\n"
            "output_new_params_fn\n"
            "output_data_fn\n"
            "num_retargets\n"
            "energy_threshold\n"
        );
        printf("\n");
        return EXIT_SUCCESS;
    }
    
    int ind=1;
    char* oldAnsatzFN = varg[ind++];
    char* oldParamsFN = varg[ind++];
    char* newAnsatzFN = varg[ind++];
    char* newParamsFN = varg[ind++];
    char* trueStateFN = varg[ind++];
    char* outParamsFN = varg[ind++];
    char* outDataFN   = varg[ind++];
    int numRetargets  = atoi(varg[ind++]);
    double energyThreshold = strtod(varg[ind++], NULL); 
    
    printf("\n");
    printf("old ansatz:\t%s\n", oldAnsatzFN);
    printf("old params:\t%s\n", oldParamsFN);
    printf("new ansatz:\t%s\n", newAnsatzFN);
    printf("new params:\t%s\n", newParamsFN);
    printf("true state:\t%s\n", trueStateFN);
    printf("out params:\t%s\n", outParamsFN);
    printf("out data:  \t%s\n", outDataFN);
    printf("num retargs:\t%d\n", numRetargets);
    printf("E threshold:\t%lf\n", energyThreshold);
    printf("\n");



   /* write simulation constants to file */

    FILE* outFile = openAssocWrite(outDataFN);
    writeDoubleToAssoc(outFile, 
        "IMAG_EVO_TIME_STEP", IMAG_EVO_TIME_STEP, OUTPUT_NUM_SIGFIGS);
    writeIntToAssoc(outFile, 
        "IMAG_EVO_NUM_ITERS", IMAG_EVO_NUM_ITERS);
    writeDoubleToAssoc(outFile, 
        "energyThreshold", energyThreshold, OUTPUT_NUM_SIGFIGS);
    writeIntToAssoc(outFile, 
        "numRetargets", numRetargets);
    writeIntToAssoc(outFile, 
        "OUTPUT_NUM_SIGFIGS", OUTPUT_NUM_SIGFIGS);
    writeStringToAssoc(outFile,
        "oldAnsatzFN", oldAnsatzFN);
    writeStringToAssoc(outFile,
        "newAnsatzFN", newAnsatzFN);
    writeStringToAssoc(outFile,
        "oldParamsFN", oldParamsFN);
    writeStringToAssoc(outFile,
        "trueStateFN", trueStateFN);
    writeStringToAssoc(outFile,
        "newParamsInitFN", newParamsFN);
    writeStringToAssoc(outFile,
        "newParamsFinalFN", outParamsFN);
    closeAssocWrite(outFile);

    
    
    /* prepare simulation structures */
        
    // simulated system info (newAnsatz is actually reversed)
    QuESTEnv questEnv = initQuESTEnv();
    Circuit oldAnsatz = loadCircuit(oldAnsatzFN);
    Circuit newAnsatz = loadCircuit(newAnsatzFN);
    double* oldParams = loadParams(oldParamsFN, oldAnsatz.numGates);
    double* newParams = loadParams(newParamsFN, newAnsatz.numGates);
    Hamiltonian hamil = loadHamiltonian(RECOMP_HAMIL_FN, questEnv);
    ParamEvolEnv paramEnv = initParamEvolEnv(
        newAnsatz.numQubits, newAnsatz.numGates, questEnv);
    
    // old ansatz initial state: |++++++1>
    QubitRegister initState = createQubitRegister(oldAnsatz.numQubits, questEnv);
    initStatePlus(initState);
    hadamard(initState, 0);
    sigmaX(initState,   0);
    
    // old ansatz final state: oldAnsatz(oldParams) |initState>
    QubitRegister finalState = createQubitRegister(oldAnsatz.numQubits, questEnv);
    cloneQubitRegister(finalState, initState);
    applyCircuit(oldAnsatz, oldParams, finalState);
    
    // targetState will be slowly adjusted to finalState (unless no retargeting)
    double targetParams[oldAnsatz.numGates];
    for (int p=0; p < oldAnsatz.numGates; p++) {
        if (numRetargets == 0)
            targetParams[p] = oldParams[p];
        else
            targetParams[p] = IMAG_EVO_INIT_TARGET_PARAM;   
    }    
    QubitRegister targetState = createQubitRegister(oldAnsatz.numQubits, questEnv);
    cloneQubitRegister(targetState, initState);
    applyCircuit(oldAnsatz, targetParams, targetState);
    
    // imaginary parameterised wavefunction
    QubitRegister paramState = createQubitRegister(newAnsatz.numQubits, questEnv);
    cloneQubitRegister(paramState, initState);
    applyCircuit(oldAnsatz, targetParams, paramState);
    applyCircuit(newAnsatz, newParams, paramState);
    
    // true state (result of true time evolution) which oldAnsatz approximates
    QubitRegister trueState = createQubitRegister(newAnsatz.numQubits, questEnv);
    loadWavefunction(trueStateFN, trueState);
    
    // result of inv(newAnsatz)|in> (to approximate oldAnatz|in> and trueState)
    QubitRegister revNewState = createQubitRegister(newAnsatz.numQubits, questEnv);
    cloneQubitRegister(revNewState, initState);
    applyInverseCircuit(newAnsatz, newParams, revNewState);
    
    
    
    /* prepare data collecting structures */
    
    // prepare output data structures
    double* fidelityWithTargetEvo = createArray(IMAG_EVO_NUM_ITERS);
    double* fidelityWithFinalEvo = createArray(IMAG_EVO_NUM_ITERS);
    double* fidelityWithTrueEvo = createArray(IMAG_EVO_NUM_ITERS);
    double* energyWithTargetEvo = createArray(IMAG_EVO_NUM_ITERS);
    double* energyWithFinalEvo = createArray(IMAG_EVO_NUM_ITERS); 
    double** newParamEvos = createNestedArray(paramEnv.numParams, IMAG_EVO_NUM_ITERS);
    double** targetParamEvos = createNestedArray(oldAnsatz.numGates, IMAG_EVO_NUM_ITERS);
    int* retargetIters = malloc(numRetargets * sizeof(int));
    for (int k=0; k < numRetargets; k++)
        retargetIters[k] = -1;
    


    /* simulate */
    
    double currEnergy = getEnergy(hamil, paramState);
    int currTarget = 0;
    
    // optimise new params so as to produce initState from targetState, via hamil
    for (int iter=0; iter < IMAG_EVO_NUM_ITERS; iter++) {
        
        // prepare the input state to imag-time, which slowly becomes the old reached state
        double distToGround = currEnergy - hamil.eigvals[0];
        if (distToGround < energyThreshold && currTarget < numRetargets) {
                        
            retargetIters[currTarget++] = iter;
            double fac = (currTarget / (double) numRetargets);
            for (int p=0; p < oldAnsatz.numGates; p++)
                targetParams[p] = fac * oldParams[p];

            cloneQubitRegister(targetState, initState);
            applyCircuit(oldAnsatz, targetParams, targetState);
            
            printf(
                "\tUpdated target params to %d/%d (%g%%)\n", 
                currTarget, numRetargets, round(fac*100));
        }
        
        // step forward in imaginary time
        evolveParamsImag(
            newParams, targetState, newAnsatz, hamil, IMAG_EVO_TIME_STEP, paramEnv);
        
        // reccord parameter evolution
        for (int p=0; p < paramEnv.numParams; p++)
            newParamEvos[p][iter] = newParams[p]; 
        for (int p=0; p < oldAnsatz.numGates; p++)
            targetParamEvos[p][iter] = targetParams[p];
            
        // measure how well newParams deconstruct targetState
        cloneQubitRegister(paramState, targetState);
        applyCircuit(newAnsatz, newParams, paramState);
        fidelityWithTargetEvo[iter] = calcFidelity(paramState, initState);
        energyWithTargetEvo[iter] = getEnergy(hamil, paramState);
                
        // measure how well newParams deconstruct finalState
        cloneQubitRegister(paramState, finalState);
        applyCircuit(newAnsatz, newParams, paramState);
        fidelityWithFinalEvo[iter] = calcFidelity(paramState, initState);
        energyWithFinalEvo[iter] = getEnergy(hamil, paramState);
        
        // measure how well newParams reconstruct trueState
        cloneQubitRegister(revNewState, initState);
        applyInverseCircuit(newAnsatz, newParams, revNewState);
        fidelityWithTrueEvo[iter] = calcFidelity(revNewState, trueState);
        
        // monitor progress
        currEnergy = energyWithTargetEvo[iter];
        if (numRetargets == 0)
            printf(
                "t=%d\tfidelity: %lf\tenergy: %lf (old) f: %g (true)\n", 
                iter, fidelityWithFinalEvo[iter], energyWithFinalEvo[iter],
                fidelityWithTrueEvo[iter]);
        else
            printf(
                "t=%d\tf: %lf  E: %lf (target) f: %lf  E: %lf (final) f: %g (true)\n",
                iter, fidelityWithTargetEvo[iter], energyWithTargetEvo[iter],
                fidelityWithFinalEvo[iter], energyWithFinalEvo[iter],
                fidelityWithTrueEvo[iter]);
    }
    
    
    
    /* write data to file */
    
    // write imag-time simulation results to file
    outFile = openAssocAppend(outDataFN);
    writeDoubleArrToAssoc(outFile, 
        "hamilSpectrum", hamil.eigvals, 2LL << (hamil.numQubits-1), 
        OUTPUT_NUM_SIGFIGS);
    writeDoubleArrToAssoc(outFile, 
        "fidelityWithTargetEvo", fidelityWithTargetEvo, 
        IMAG_EVO_NUM_ITERS, OUTPUT_NUM_SIGFIGS);
    writeDoubleArrToAssoc(outFile, 
        "fidelityWithFinalEvo", fidelityWithFinalEvo, 
        IMAG_EVO_NUM_ITERS, OUTPUT_NUM_SIGFIGS);
    writeDoubleArrToAssoc(outFile, 
        "fidelityWithTrueEvo", fidelityWithTrueEvo, 
        IMAG_EVO_NUM_ITERS, OUTPUT_NUM_SIGFIGS);
    writeDoubleArrToAssoc(outFile, 
        "energyWithTargetEvo", energyWithTargetEvo, 
        IMAG_EVO_NUM_ITERS, OUTPUT_NUM_SIGFIGS);
    writeDoubleArrToAssoc(outFile, 
        "energyWithFinalEvo", energyWithFinalEvo, 
        IMAG_EVO_NUM_ITERS, OUTPUT_NUM_SIGFIGS);
    writeDoubleArrToAssoc(outFile,
        "oldParams", oldParams, oldAnsatz.numGates, 
        OUTPUT_NUM_SIGFIGS);
    writeNestedDoubleListToAssoc(outFile, 
        "newParamEvos", newParamEvos, 
        2, (int []) {newAnsatz.numGates, IMAG_EVO_NUM_ITERS}, 
        OUTPUT_NUM_SIGFIGS);
    writeNestedDoubleListToAssoc(outFile, 
        "targetParamEvos", targetParamEvos, 
        2, (int []) {oldAnsatz.numGates, IMAG_EVO_NUM_ITERS}, 
        OUTPUT_NUM_SIGFIGS);
    writeIntArrToAssoc(outFile, 
        "retargetIters", retargetIters, currTarget);
    closeAssocAppend(outFile);
    
    // record final new param values seperately
    char header[500];
    sprintf(header,
        "# Params after %d iters (dt=%lf) in imaginary-time of %s under %s, "
        "recompiling %s with parameters %s",
        IMAG_EVO_NUM_ITERS, IMAG_EVO_TIME_STEP, newAnsatzFN, RECOMP_HAMIL_FN, 
        oldAnsatzFN, oldParamsFN);
    writeParams(outParamsFN, newParams, newAnsatz.numGates, header);



    /* free structures */
    
    // data collecting
    free(fidelityWithTargetEvo);
    free(fidelityWithFinalEvo);
    free(fidelityWithTrueEvo);
    free(energyWithTargetEvo);
    free(energyWithFinalEvo);
    free(retargetIters);
    freeNestedArray(newParamEvos, paramEnv.numParams);
    freeNestedArray(targetParamEvos, oldAnsatz.numGates);
    
    // simulation
    free(oldParams);
    free(newParams);
    freeCircuit(oldAnsatz);
    freeCircuit(newAnsatz);
    freeHamiltonian(hamil, questEnv);
    destroyQubitRegister(paramState, questEnv);
    destroyQubitRegister(targetState, questEnv);
    destroyQubitRegister(finalState, questEnv);
    destroyQubitRegister(initState, questEnv);
    destroyQubitRegister(revNewState, questEnv);
    destroyQubitRegister(trueState, questEnv);
    closeParamEvolEnv(paramEnv, questEnv);
    closeQuESTEnv(questEnv);    
    return EXIT_SUCCESS;
}