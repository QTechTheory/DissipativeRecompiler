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
#define IMAG_SOLVER_METHOD 0        // 0 for TSVD, 1 for Tikhonov

#define ELIM_PARAM_MAX 1E-6
#define PARAM_CHANGE_MAX 1E-1
#define ENERGY_INCREASE_TOL_FAC 2

#define MAX_TOTAL_ITERS 5000
#define OUTPUT_NUM_SIGFIGS 10



int getAbsMinInd(double* list, int len, int* include) {
    
    int ind = 0;
    
    for (int i=1; i < len; i++)
        if (fabs(list[i]) < fabs(list[ind]) && include[i])
            ind = i;
    
    return ind;
}



int main(int narg, char *varg[]) {
    
    if (narg != 1 + 7) {
        printf("\nargs:\n"
            "old_ansatz_fn\n"
            "old_params_fn\n"
            "new_ansatz_fn\n"
            "new_params_fn\n"
            "true_state_fn\n"
            "out_params_fn\n"
            "output_fn\n\n");
        exit(EXIT_SUCCESS);
    }
    
    
    
    int ind=1;
    char* oldAnsatzFN = varg[ind++];
    char* oldParamsFN = varg[ind++];
    char* newAnsatzFN = varg[ind++];
    char* newParamsFN = varg[ind++];
    char* trueStateFN = varg[ind++];
    char* outParamsFN = varg[ind++];
    char* outDataFN   = varg[ind++];
    
    printf("\n");
    printf("old ansatz:\t%s\n", oldAnsatzFN);
    printf("old params:\t%s\n", oldParamsFN);
    printf("new ansatz:\t%s\n", newAnsatzFN);
    printf("new params:\t%s\n", newParamsFN);
    printf("true state:\t%s\n", trueStateFN);
    printf("out params:\t%s\n", outParamsFN);
    printf("out data:  \t%s\n", outDataFN);
    
    
    
    /* write simulation constants to file */
    
    FILE* outFile = openAssocWrite(outDataFN);
    writeStringToAssoc(outFile,
        "RECOMP_HAMIL_FN", RECOMP_HAMIL_FN);
    writeDoubleToAssoc(outFile, 
        "IMAG_EVO_TIME_STEP", IMAG_EVO_TIME_STEP, OUTPUT_NUM_SIGFIGS);
    writeIntToAssoc(outFile,
        "IMAG_SOLVER_METHOD", IMAG_SOLVER_METHOD);
    writeIntToAssoc(outFile, 
        "OUTPUT_NUM_SIGFIGS", OUTPUT_NUM_SIGFIGS);
    writeIntToAssoc(outFile, 
        "MAX_TOTAL_ITERS", MAX_TOTAL_ITERS);
    writeDoubleToAssoc(outFile, 
        "ELIM_PARAM_MAX", ELIM_PARAM_MAX, OUTPUT_NUM_SIGFIGS);
    writeDoubleToAssoc(outFile, 
        "PARAM_CHANGE_MAX", PARAM_CHANGE_MAX, OUTPUT_NUM_SIGFIGS);
    writeDoubleToAssoc(outFile, 
        "ENERGY_INCREASE_TOL_FAC", ENERGY_INCREASE_TOL_FAC, OUTPUT_NUM_SIGFIGS);
    writeStringToAssoc(outFile,
        "oldAnsatzFN", oldAnsatzFN);
    writeStringToAssoc(outFile,
        "newAnsatzFN", newAnsatzFN);
    writeStringToAssoc(outFile,
        "trueStateFN", trueStateFN);
    writeStringToAssoc(outFile,
        "oldParamsFN", oldParamsFN);
    writeStringToAssoc(outFile,
        "newParamsInitFN", newParamsFN);
    writeStringToAssoc(outFile,
        "newParamsFinalFN", outParamsFN);    
    closeAssocWrite(outFile);


    
    /* prepare simulation structures */
    
    QuESTEnv questEnv = initQuESTEnv();
    Circuit oldAnsatz = loadCircuit(oldAnsatzFN);
    Circuit newAnsatz = loadCircuit(newAnsatzFN);
    double* oldParams = loadParams(oldParamsFN, oldAnsatz.numGates);
    double* newParams = loadParams(newParamsFN, newAnsatz.numGates);
    Hamiltonian hamil = loadHamiltonian(RECOMP_HAMIL_FN, questEnv);
    ParamEvolEnv paramEnv = initParamEvolEnv(newAnsatz.numQubits, newAnsatz.numGates, questEnv);    
    
    // initial state: |++++++1> (ground state of hamil)
    QubitRegister initState = createQubitRegister(newAnsatz.numQubits, questEnv);
    initStatePlus(initState);
    hadamard(initState, 0);
    sigmaX(initState,   0);
    
    // state to be 'uncompiled': oldAnsatz(oldParams)|initState>
    QubitRegister oldState = createQubitRegister(oldAnsatz.numQubits, questEnv);
    cloneQubitRegister(oldState, initState);
    applyCircuit(oldAnsatz, oldParams, oldState);
    
    // paramState = newAnsatz(newParams) oldAnsatz(oldParams) |initState>
    QubitRegister paramState = createQubitRegister(newAnsatz.numQubits, questEnv);
    cloneQubitRegister(paramState, oldState);
    applyCircuit(newAnsatz, newParams, paramState);
    
    // backup for restoring previously low-energy parameters
    double* newParamsBackup = createArray(paramEnv.numParams);
    cloneArray(newParamsBackup, newParams, paramEnv.numParams);
    
    // 1 means param active, 0 means it's been shut off
    int* newParamStatus = malloc(paramEnv.numParams * sizeof *newParamStatus);
    for (int i=0; i < paramEnv.numParams; i++)
        newParamStatus[i] = 1;
    
    // true state (result of true time evolution)
    QubitRegister trueState = createQubitRegister(newAnsatz.numQubits, questEnv);
    loadWavefunction(trueStateFN, trueState);
    
    // result of inv(newAnsatz)|in> (to approximate oldAnatz|in> and trueState)
    QubitRegister revNewState = createQubitRegister(newAnsatz.numQubits, questEnv);
    cloneQubitRegister(revNewState, initState);
    applyInverseCircuit(newAnsatz, newParams, revNewState);
    
    
    /* prepare data recording structures */
    int totalIters = 0;
    int totalParamChanges = 0; // excludes final failed
    int* itersOfParamChanges = malloc(paramEnv.numParams * sizeof* itersOfParamChanges);
    int* indsOfParamChanges = malloc(paramEnv.numParams * sizeof* indsOfParamChanges);
    double* fidelityInitEvo = createArray(MAX_TOTAL_ITERS);
    double* fidelityTrueEvo = createArray(MAX_TOTAL_ITERS);
    double* energyEvo = createArray(MAX_TOTAL_ITERS);
    double** paramsEvo = createNestedArray(paramEnv.numParams, MAX_TOTAL_ITERS);
    
    
    
    /* perform gate elimination */
    
    // the starting energy determines stopping criteria 
    double currEnergy = getEnergy(hamil, paramState);
    double currDist = currEnergy - hamil.eigvals[0];
    double origDist = currDist;
    double currInitFidelity = calcFidelity(paramState, initState);
    double currTrueFidelity = calcFidelity(revNewState, trueState);
    
    // record initial energy / params
    energyEvo[0] = currEnergy;
    fidelityInitEvo[0] = currInitFidelity;
    fidelityTrueEvo[0] = currTrueFidelity;
    for (int p=0; p < paramEnv.numParams; p++)
        paramsEvo[p][0] = newParams[p];
    
    printf("Initial energy: %lf (%lf from ground)\n", currEnergy, origDist);
    
    // pick the smallest parameter gate
    int paramInd = getAbsMinInd(newParams, paramEnv.numParams, newParamStatus);
    itersOfParamChanges[0] = 0;
    indsOfParamChanges[0] = paramInd;
    
    // keep eliminating until all gates sufficiently small, or energy grew too high
    while (
        fabs(newParams[paramInd]) > ELIM_PARAM_MAX &&
        currDist < origDist * ENERGY_INCREASE_TOL_FAC
    ) {
        
        printf("Attempting to eliminate gate %d (val=%lf)\n", paramInd, newParams[paramInd]);
        
        // remember param config in case we need to restore it
        cloneArray(newParamsBackup, newParams, paramEnv.numParams);
    
        // continue constrained imag-evol until param is zero or energy grows badly
        while (
            fabs(newParams[paramInd]) > ELIM_PARAM_MAX && 
            currDist < origDist * ENERGY_INCREASE_TOL_FAC
        ) {
            // safety check
            if (totalIters >= MAX_TOTAL_ITERS-1) {
                printf("ERROR! Reached MAX_TOTAL_ITERS=%d! Enlargen this to allocate more data space",
                    MAX_TOTAL_ITERS);
                exit(EXIT_FAILURE);
            }

            // populate the matrices M dnewParams = V
            populateMatrices(paramEnv, oldState, newAnsatz, newParams, hamil, 0);
        
            // modify the matrices to push param toward 0
            double forcedChange = -newParams[paramInd] / IMAG_EVO_TIME_STEP;
            if (fabs(forcedChange) > PARAM_CHANGE_MAX)
                forcedChange = (forcedChange > 0)? PARAM_CHANGE_MAX : -PARAM_CHANGE_MAX;

            // by moving the the paramInd column (times change) into RHS
            for (int i=0; i < paramEnv.numParams; i++) {
                double oldVecVal = gsl_vector_get(paramEnv.varVector, i);
                double oldMatVal = gsl_matrix_get(paramEnv.varMatrix, i, paramInd);
                double newVecVal = oldVecVal - oldMatVal * forcedChange;
                double newMatVal = 0;
                gsl_vector_set(paramEnv.varVector, i, newVecVal);
                gsl_matrix_set(paramEnv.varMatrix, i, paramInd, newMatVal);
            }
            
            // also constrain that all disabled params don't change
            for (int p=0; p < paramEnv.numParams; p++) {
                
                // by zeroing their columns
                if (newParamStatus[p] == 0)
                    for (int i=0; i < paramEnv.numParams; i++)
                        gsl_matrix_set(paramEnv.varMatrix, i, p, 0);
            }
            
            // solve the modified equations
            if (IMAG_SOLVER_METHOD == 0)
                solveViaTSVD(paramEnv);
            else
                solveViaTikhonov(paramEnv);
            
            // update only the enabled parameters
            for (int p=0; p < paramEnv.numParams; p++) {
                if (p == paramInd)
                    newParams[p] += forcedChange * IMAG_EVO_TIME_STEP;
                else if (newParamStatus[p])
                    newParams[p] += IMAG_EVO_TIME_STEP * gsl_vector_get(paramEnv.paramChange, p);
            }   
            
            // create the new parameterised state
            cloneQubitRegister(paramState, oldState);
            applyCircuit(newAnsatz, newParams, paramState);
            
            // (and the state when forward applied)
            cloneQubitRegister(revNewState, initState);
            applyInverseCircuit(newAnsatz, newParams, revNewState);
                        
            // recheck energy
            currEnergy = getEnergy(hamil, paramState);
            currDist = currEnergy - hamil.eigvals[0];
            currInitFidelity = calcFidelity(paramState, initState);
            currTrueFidelity = calcFidelity(revNewState, trueState);
            
            // record data
            totalIters++;
            energyEvo[totalIters] = currEnergy;
            fidelityInitEvo[totalIters] = currInitFidelity;
            fidelityTrueEvo[totalIters] = currTrueFidelity;
            for (int p=0; p < paramEnv.numParams; p++)
                paramsEvo[p][totalIters] = newParams[p];
            
            
            printf(
                "param[%d]=%g,\tdist=%g,\tfid|in>=%g,\tfid|true>=%g\n", 
                paramInd, newParams[paramInd], currDist, currInitFidelity, currTrueFidelity);
        }
        
        // if energy rose too much, restore the original parameter
        if (currDist >= origDist * ENERGY_INCREASE_TOL_FAC) {
            
            printf(
                "Oh no! Energy gap became %lf > %d * %lf = %lf\n", 
                currDist, ENERGY_INCREASE_TOL_FAC, origDist, 
                origDist * ENERGY_INCREASE_TOL_FAC
            );
                
            cloneArray(newParams, newParamsBackup, paramEnv.numParams);
            cloneQubitRegister(paramState, oldState);
            applyCircuit(newAnsatz, newParams, paramState);
            currEnergy = getEnergy(hamil, paramState);
            currDist = currEnergy - hamil.eigvals[0];
            
            break;
            
            // @TODO IF YOU ALLOWED ITERATION TO CONTINUE HERE, MARK THIS ITERATION HAD A RESET IN OUTPUT
        }
        
        // otherwise record success and zero more parameters!
        newParamStatus[paramInd] = 0;

        printf(
            "Yay! param[%d] = %lf < %lf\n", 
            paramInd, newParams[paramInd], ELIM_PARAM_MAX);
            
        // pick the next smallest param to remove
        paramInd = getAbsMinInd(newParams, paramEnv.numParams, newParamStatus);
        totalParamChanges++;
        itersOfParamChanges[totalParamChanges] = totalIters;
        indsOfParamChanges[totalParamChanges] = paramInd;
    }
    
    
    
    /* write data to file */
    
    // simulation data
    outFile = openAssocAppend(outDataFN);
    writeDoubleArrToAssoc(outFile, 
        "hamilSpectrum", hamil.eigvals, 
        2LL << (hamil.numQubits-1), OUTPUT_NUM_SIGFIGS);
    writeIntToAssoc(outFile, 
        "totalIters", totalIters);
    writeIntToAssoc(outFile, 
        "totalParamChanges", totalParamChanges);
    writeDoubleToAssoc(outFile,
        "origDist", origDist, OUTPUT_NUM_SIGFIGS);
    writeIntArrToAssoc(outFile, 
        "itersOfParamChanges", itersOfParamChanges, totalParamChanges);
    writeIntArrToAssoc(outFile, 
        "indsOfParamChanges", indsOfParamChanges, totalParamChanges);
    writeDoubleArrToAssoc(outFile, 
        "fidelityInitEvo", fidelityInitEvo, totalIters, OUTPUT_NUM_SIGFIGS);
    writeDoubleArrToAssoc(outFile, 
        "fidelityTrueEvo", fidelityTrueEvo, totalIters, OUTPUT_NUM_SIGFIGS);
    writeDoubleArrToAssoc(outFile, 
        "energyEvo", energyEvo, totalIters, OUTPUT_NUM_SIGFIGS);
    writeNestedDoubleListToAssoc(outFile, 
        "paramsEvo", paramsEvo, 2, (int []) {paramEnv.numParams, totalIters}, 
        OUTPUT_NUM_SIGFIGS);
    closeAssocAppend(outFile);

    // record final new param values seperately
    char header[500];
    sprintf(header,
        "# Params after elimination via %d iters (dt=%lf) in imaginary-time of %s under %s, "
        "starting from %s and eliminating %d.",
        totalIters, IMAG_EVO_TIME_STEP, newAnsatzFN, RECOMP_HAMIL_FN, 
        newParamsFN, totalParamChanges-1);
    writeParams(outParamsFN, newParams, newAnsatz.numGates, header);



    /* free structures */
    
    // data collecting
    free(itersOfParamChanges);
    free(indsOfParamChanges);
    free(energyEvo);
    free(fidelityInitEvo);
    free(fidelityTrueEvo);
    freeNestedArray(paramsEvo, paramEnv.numParams);

    // simulation
    free(oldParams);
    free(newParams);
    free(newParamStatus);
    freeCircuit(oldAnsatz);
    freeCircuit(newAnsatz);
    freeHamiltonian(hamil, questEnv);
    destroyQubitRegister(initState, questEnv);
    destroyQubitRegister(oldState, questEnv);
    destroyQubitRegister(paramState, questEnv);
    destroyQubitRegister(revNewState, questEnv);
    destroyQubitRegister(trueState, questEnv);
    closeParamEvolEnv(paramEnv, questEnv);
    closeQuESTEnv(questEnv);    
    return EXIT_SUCCESS;
}