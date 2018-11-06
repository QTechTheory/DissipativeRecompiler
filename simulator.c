#include <QuEST.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

#include "circuitloader.h"
#include "hamiltonianloader.h"
#include "paramevolver.h"
#include "trotterevolver.h"
#include "trueevolver.h"
#include "utilities.h"
#include "mmaformatter.h"



#define TROT_CIRC_FN "data/trotter_cycle.txt"
#define SIM_HAMIL_FN "data/hamil_sim.txt"

#define REAL_EVO_TIME_STEP 0.0025
#define REAL_EVO_NUM_ITERS 700      // fid=0.995 at iter~700
#define REAL_EVO_TROT_MODE 1        // 0=fixed, 1=alt.rev. 2=random
#define REAL_EVO_TROT_CYCLES 0      // 0=no trotter simulation

#define RAND_TROT_SEED 123
#define OUTPUT_NUM_SIGFIGS 10



int main(int narg, char *varg[]) {

    /* collect arguments */
    
    if (narg != 1 + 5 && narg != 1 + 6) {
        printf("\nargs:\n"
            "ansatz_fn\n"
            "init_param_fn\n"
            "output_param_fn\n"
            "output_true_wavef_fn\n"
            "output_data_fn\n"
            "[input_true_init_wavef_fn]\n");
        printf("\n");
        exit(EXIT_SUCCESS);
    }
    
    int ind=1;
    char* inAnsatzFN = varg[ind++];
    char* inParamFN  = varg[ind++];
    char* outParamFN = varg[ind++];
    char* outWavefFN = varg[ind++];
    char* outDataFN  = varg[ind++];
    char* inWavefFN  = (narg == 1 + 6)? varg[ind++] : NULL;
    
    printf("\n");
    if (inWavefFN != NULL)
        printf("in wavef:\t%s\n", inWavefFN);
    printf("in params:\t%s\n",inParamFN);
    printf("in ansatz: \t%s\n", inAnsatzFN);
    printf("out params:\t%s\n", outParamFN);
    printf("out wavef: \t%s\n", outWavefFN);
    printf("out data:\t%s\n\n", outDataFN);
    


    /* write simulation constants to file */
    
    FILE* outFile = openAssocWrite(outDataFN);
    writeDoubleToAssoc(outFile, 
        "REAL_EVO_TIME_STEP", REAL_EVO_TIME_STEP, OUTPUT_NUM_SIGFIGS);
    writeIntToAssoc(outFile, 
        "REAL_EVO_NUM_ITERS", REAL_EVO_NUM_ITERS);
    writeIntToAssoc(outFile, 
        "REAL_EVO_TROT_MODE", REAL_EVO_TROT_MODE);
    writeIntToAssoc(outFile, 
        "REAL_EVO_TROT_CYCLES", REAL_EVO_TROT_CYCLES);
    writeIntToAssoc(outFile, 
        "RAND_TROT_SEED", RAND_TROT_SEED);
    writeIntToAssoc(outFile, 
        "OUTPUT_NUM_SIGFIGS", OUTPUT_NUM_SIGFIGS);
    writeStringToAssoc(outFile,
        "inAnsatzFN", inAnsatzFN);
    writeStringToAssoc(outFile,
        "inParamFN", inParamFN);
    writeStringToAssoc(outFile,
        "outParamFN", outParamFN);
    writeStringToAssoc(outFile,
        "outWavefFN", outWavefFN);
    if (inWavefFN != NULL)
        writeStringToAssoc(outFile,
            "inTrueInitWavefFN", inWavefFN);
    closeAssocWrite(outFile);
    
    
    
    /* prepare simulation structures */
    
    // seed random Trotter ordering
    srand(RAND_TROT_SEED);
    
    // simulated system info
    QuESTEnv questEnv = initQuESTEnv();
    Circuit ansatz = loadCircuit(inAnsatzFN);
    Hamiltonian hamil = loadHamiltonian(SIM_HAMIL_FN, questEnv);
    TrueEvolEnv trueEnv = initTrueEvolEnv(hamil, REAL_EVO_TIME_STEP, questEnv);
    ParamEvolEnv paramEnv = initParamEvolEnv(ansatz.numQubits, ansatz.numGates, questEnv);    
    
    // initial state |++++++1>
    QubitRegister initState = createQubitRegister(ansatz.numQubits, questEnv);
    initStatePlus(initState);
    hadamard(initState, 0);
    sigmaX(initState,   0);
    
    // true evolution wavefunction
    QubitRegister trueState = createQubitRegister(ansatz.numQubits, questEnv);
    if (inWavefFN == NULL)
        cloneQubitRegister(trueState, initState);
    else
        loadWavefunction(inWavefFN, trueState);
    
    // trotter wavefunction
    QubitRegister trotterState = createQubitRegister(ansatz.numQubits, questEnv);
    Circuit trotterCirc = loadCircuit(TROT_CIRC_FN);
    double trotterAngles[trotterCirc.numGates];
    
    // parameterised wavefunction
    QubitRegister paramState = createQubitRegister(ansatz.numQubits, questEnv);
    double* params = loadParams(inParamFN, ansatz.numGates);



    /* prepare data collecting structures */

    double* variationalFidelityEvo = createArray(REAL_EVO_NUM_ITERS);
    double* trotterFidelityEvo = createArray(REAL_EVO_NUM_ITERS);
    double** paramEvos = createNestedArray(paramEnv.numParams, REAL_EVO_NUM_ITERS);
    double** trotterAngleEvos = createNestedArray(trotterCirc.numGates, REAL_EVO_NUM_ITERS);
    
    double* residualEvo = createArray(REAL_EVO_NUM_ITERS);
    double* energyEvo = createArray(REAL_EVO_NUM_ITERS);
    
    
    /* simulate */
    
    for (int iter=0; iter < REAL_EVO_NUM_ITERS; iter++) {

        // get true wavef
        evolveTrueState(trueState, hamil, trueEnv);
    
        // get trotter wavef
        if (REAL_EVO_TROT_CYCLES > 0) {
            double trotterTime = (iter+1) * REAL_EVO_TIME_STEP;
            produceTrotterState(
                initState, trotterState, hamil, trotterCirc, trotterAngles, trotterTime, 
                REAL_EVO_TROT_CYCLES, REAL_EVO_TROT_MODE);    
     
            trotterFidelityEvo[iter] = calcFidelity(trueState, trotterState);
            for (int p=0; p < trotterCirc.numGates; p++)
                trotterAngleEvos[p][iter] = trotterAngles[p];
        }
    
        // get parameterised wavef
        residualEvo[iter] = evolveParamsReal(params, initState, ansatz, hamil, REAL_EVO_TIME_STEP, paramEnv);
        cloneQubitRegister(paramState, initState);
        applyCircuit(ansatz, params, paramState);
        
        energyEvo[iter] = getEnergy(hamil, paramState);
        variationalFidelityEvo[iter] = calcFidelity(trueState, paramState);
        for (int p=0; p < paramEnv.numParams; p++)
            paramEvos[p][iter] = params[p];
    
        // monitor progress
        if (REAL_EVO_TROT_CYCLES > 0)
            printf(
                "t=%d\ttrotter fidelity: %lf\tvariational fidelity: %lf\n", iter, 
                trotterFidelityEvo[iter], variationalFidelityEvo[iter]);
        else
            printf(
                "t=%d\tvariational fidelity: %lf, resid: %g,\tenergy: %lf\n", 
                iter, variationalFidelityEvo[iter], residualEvo[iter], energyEvo[iter]);     
    }
    
    
    
    /* write data to file */
    
    outFile = openAssocAppend(outDataFN);
    writeDoubleArrToAssoc(outFile, 
        "variationalFidelityEvo", variationalFidelityEvo,
        REAL_EVO_NUM_ITERS, OUTPUT_NUM_SIGFIGS);
    writeDoubleArrToAssoc(outFile, 
        "trotterFidelityEvo", trotterFidelityEvo,
        REAL_EVO_NUM_ITERS, OUTPUT_NUM_SIGFIGS);
    writeDoubleArrToAssoc(outFile,
        "energyEvo", energyEvo,
        REAL_EVO_NUM_ITERS, OUTPUT_NUM_SIGFIGS);
    writeDoubleArrToAssoc(outFile,
        "residualEvo", residualEvo,
        REAL_EVO_NUM_ITERS, OUTPUT_NUM_SIGFIGS);
    writeNestedDoubleListToAssoc(outFile, 
        "paramEvos", paramEvos, 
        2, (int []) {paramEnv.numParams, REAL_EVO_NUM_ITERS}, OUTPUT_NUM_SIGFIGS);
    writeNestedDoubleListToAssoc(outFile, 
        "trotterAngleEvos", trotterAngleEvos, 
        2, (int []) {trotterCirc.numGates, REAL_EVO_NUM_ITERS}, OUTPUT_NUM_SIGFIGS);
    writeDoubleArrToAssoc(outFile, 
        "hamilSpectrum", hamil.eigvals,
        2LL << (hamil.numQubits-1), OUTPUT_NUM_SIGFIGS);
    closeAssocWrite(outFile);
    
    // record final param values separately
    char header[500];
    sprintf(header,
        "# Params after %d iters (dt=%lf) in real-time of %s under %s, starting from %s", 
        REAL_EVO_NUM_ITERS, REAL_EVO_TIME_STEP, inAnsatzFN, SIM_HAMIL_FN, inParamFN);
    writeParams(outParamFN, params, ansatz.numGates, header);
    
    // record final true wavefunction seperately
    sprintf(header,
        "# True wavefunction after %d iters (dt=%lf) in real-time under %s", 
        REAL_EVO_NUM_ITERS, REAL_EVO_TIME_STEP, SIM_HAMIL_FN);
    writeWavefunction(outWavefFN, trueState, header);
    
    
    /* free structures */
    
    // data collecting 
    free(residualEvo);
    free(energyEvo);
    free(variationalFidelityEvo);
    free(trotterFidelityEvo);
    freeNestedArray(paramEvos, paramEnv.numParams);
    freeNestedArray(trotterAngleEvos, trotterCirc.numGates);

    // simulation
    free(params);
    freeCircuit(ansatz);
    freeCircuit(trotterCirc);
    freeHamiltonian(hamil, questEnv);
    destroyQubitRegister(trueState, questEnv);
    destroyQubitRegister(initState, questEnv);
    destroyQubitRegister(paramState, questEnv);
    destroyQubitRegister(trotterState, questEnv);
    closeParamEvolEnv(paramEnv, questEnv);
    closeTrueEvolEnv(trueEnv, questEnv);
    closeQuESTEnv(questEnv);
    return EXIT_SUCCESS;
}