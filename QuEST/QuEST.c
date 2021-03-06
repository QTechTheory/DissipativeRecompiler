// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Implements the QuEST.h API (and some debugging functions) in a hardware-agnostic way, 
 * for both pure and mixed states. These functions mostly wrap hardware-specific functions,
 * and should never call eachother.
 *
 * Density matrices rho of N qubits are flattened to appear as state-vectors |s> of 2N qubits.
 * Operations U rho U^dag are implemented as U^* U |s> and make use of the pure state backend,
 * and often don't need to explicitly compute U^*.
 */

// @TODO: 4 noise functions on CPU

# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_internal.h"
# include "QuEST_validation.h"
# include "QuEST_ops.h"
# include "QuEST_qasm.h"

#ifdef __cplusplus
extern "C" {
#endif

QubitRegister createQubitRegister(int numQubits, QuESTEnv env) {
    validateCreateNumQubits(numQubits, __func__);
    
    QubitRegister qureg;
    statevec_createQubitRegister(&qureg, numQubits, env);
    qureg.isDensityMatrix = 0;
    qureg.numQubitsRepresented = numQubits;
    qureg.numQubitsInStateVec = numQubits;
    
    qasm_setup(&qureg);
    initStateZero(qureg);
    return qureg;
}

QubitRegister createDensityQubitRegister(int numQubits, QuESTEnv env) {
    validateCreateNumQubits(numQubits, __func__);
    
    QubitRegister qureg;
    statevec_createQubitRegister(&qureg, 2*numQubits, env);
    qureg.isDensityMatrix = 1;
    qureg.numQubitsRepresented = numQubits;
    qureg.numQubitsInStateVec = 2*numQubits;
    
    qasm_setup(&qureg);
    initStateZero(qureg);
    return qureg;
}

void destroyQubitRegister(QubitRegister qureg, QuESTEnv env) {
    statevec_destroyQubitRegister(qureg, env);
    qasm_free(qureg);
}

void startRecordingQASM(QubitRegister qureg) {
    qasm_startRecording(qureg);
}

void stopRecordingQASM(QubitRegister qureg) {
    qasm_stopRecording(qureg);
}

void clearRecordedQASM(QubitRegister qureg) {
    qasm_clearRecorded(qureg);
}

void printRecordedQASM(QubitRegister qureg) {
    qasm_printRecorded(qureg);
}

void writeRecordedQASMToFile(QubitRegister qureg, char* filename) {
    int success = qasm_writeRecordedToFile(qureg, filename);
    validateFileOpened(success, __func__);
}

void initStateZero(QubitRegister qureg) {
    statevec_initStateZero(qureg); // valid for both statevec and density matrices
    
    qasm_recordInitZero(qureg);
}

void initStatePlus(QubitRegister qureg) {
    if (qureg.isDensityMatrix)
        densmatr_initStatePlus(qureg);
    else
        statevec_initStatePlus(qureg);
    
    qasm_recordInitPlus(qureg);
}

void initClassicalState(QubitRegister qureg, long long int stateInd) {
    validateStateIndex(qureg, stateInd, __func__);
    
    if (qureg.isDensityMatrix)
        densmatr_initClassicalState(qureg, stateInd);
    else
        statevec_initClassicalState(qureg, stateInd);
    
    qasm_recordInitClassical(qureg, stateInd);
    
}

void initStateFromAmps(QubitRegister qureg, long long int startInd, REAL* reals, REAL* imags, long long int numAmps) {
    validateStateVecQureg(qureg, __func__);
    validateNumAmps(qureg, startInd, numAmps, __func__);
    
    statevec_initStateFromAmps(qureg, startInd, reals, imags, numAmps);
    
    qasm_recordComment(qureg, "(Initialising state from amplitude arrays not encoded)");
}

void hadamard(QubitRegister qureg, const int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_hadamard(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_hadamard(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_HADAMARD, targetQubit);
}

void rotateX(QubitRegister qureg, const int targetQubit, REAL angle) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_rotateX(qureg, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        statevec_rotateX(qureg, targetQubit+qureg.numQubitsRepresented, -angle);
    }
    
    qasm_recordParamGate(qureg, GATE_ROTATE_X, targetQubit, angle);
}

void rotateY(QubitRegister qureg, const int targetQubit, REAL angle) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_rotateY(qureg, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        statevec_rotateY(qureg, targetQubit+qureg.numQubitsRepresented, angle);
    }
    
    qasm_recordParamGate(qureg, GATE_ROTATE_Y, targetQubit, angle);
}

void rotateZ(QubitRegister qureg, const int targetQubit, REAL angle) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_rotateZ(qureg, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        statevec_rotateZ(qureg, targetQubit+qureg.numQubitsRepresented, -angle);
    }
    
    qasm_recordParamGate(qureg, GATE_ROTATE_Z, targetQubit, angle);
}

void controlledRotateX(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledRotateX(qureg, controlQubit, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledRotateX(qureg, controlQubit+shift, targetQubit+shift, -angle);
    }
    
    qasm_recordControlledParamGate(qureg, GATE_ROTATE_X, controlQubit, targetQubit, angle);
}

void controlledRotateY(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledRotateY(qureg, controlQubit, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledRotateY(qureg, controlQubit+shift, targetQubit+shift, angle); // rotateY is real
    }

    qasm_recordControlledParamGate(qureg, GATE_ROTATE_Y, controlQubit, targetQubit, angle);
}

void controlledRotateZ(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledRotateZ(qureg, controlQubit, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledRotateZ(qureg, controlQubit+shift, targetQubit+shift, -angle);
    }
    
    qasm_recordControlledParamGate(qureg, GATE_ROTATE_Z, controlQubit, targetQubit, angle);
}

void unitary(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u) {
    validateTarget(qureg, targetQubit, __func__);
    validateUnitaryMatrix(u, __func__);
    
    statevec_unitary(qureg, targetQubit, u);
    if (qureg.isDensityMatrix) {
        statevec_unitary(qureg, targetQubit+qureg.numQubitsRepresented, getConjugateMatrix(u));
    }
    
    qasm_recordUnitary(qureg, u, targetQubit);
}

void controlledUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    validateUnitaryMatrix(u, __func__);
    
    statevec_controlledUnitary(qureg, controlQubit, targetQubit, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledUnitary(qureg, controlQubit+shift, targetQubit+shift, getConjugateMatrix(u));
    }
    
    qasm_recordControlledUnitary(qureg, u, controlQubit, targetQubit);
}

void multiControlledUnitary(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u) {
    validateMultiControlsTarget(qureg, controlQubits, numControlQubits, targetQubit, __func__);
    validateUnitaryMatrix(u, __func__);
    
    statevec_multiControlledUnitary(qureg, controlQubits, numControlQubits, targetQubit, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(controlQubits, numControlQubits, shift);
        statevec_multiControlledUnitary(qureg, controlQubits, numControlQubits, targetQubit+shift, getConjugateMatrix(u));
        shiftIndices(controlQubits, numControlQubits, -shift);
    }
    
    qasm_recordMultiControlledUnitary(qureg, u, controlQubits, numControlQubits, targetQubit);
}

void compactUnitary(QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta) {
    validateTarget(qureg, targetQubit, __func__);
    validateUnitaryComplexPair(alpha, beta, __func__);
    
    statevec_compactUnitary(qureg, targetQubit, alpha, beta);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_compactUnitary(qureg, targetQubit+shift, getConjugateScalar(alpha), getConjugateScalar(beta));
    }

    qasm_recordCompactUnitary(qureg, alpha, beta, targetQubit);
}

void controlledCompactUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    validateUnitaryComplexPair(alpha, beta, __func__);
    
    statevec_controlledCompactUnitary(qureg, controlQubit, targetQubit, alpha, beta);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledCompactUnitary(qureg, 
            controlQubit+shift, targetQubit+shift, 
            getConjugateScalar(alpha), getConjugateScalar(beta));
    }
    
    qasm_recordControlledCompactUnitary(qureg, alpha, beta, controlQubit, targetQubit);
}

void sigmaX(QubitRegister qureg, const int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_sigmaX(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_sigmaX(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_SIGMA_X, targetQubit);
}

void sigmaY(QubitRegister qureg, const int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_sigmaY(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_sigmaYConj(qureg, targetQubit + qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_SIGMA_Y, targetQubit);
}

void sigmaZ(QubitRegister qureg, const int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_sigmaZ(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_sigmaZ(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_SIGMA_Z, targetQubit);
}

void sGate(QubitRegister qureg, const int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_sGate(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_sGateConj(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_S, targetQubit);
}

void tGate(QubitRegister qureg, const int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_tGate(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_tGateConj(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_T, targetQubit);
}

void phaseShift(QubitRegister qureg, const int targetQubit, REAL angle) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_phaseShift(qureg, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        statevec_phaseShift(qureg, targetQubit+qureg.numQubitsRepresented, -angle);
    }
    
    qasm_recordParamGate(qureg, GATE_PHASE_SHIFT, targetQubit, angle);
}

void controlledPhaseShift(QubitRegister qureg, const int idQubit1, const int idQubit2, REAL angle) {
    validateControlTarget(qureg, idQubit1, idQubit2, __func__);
    
    statevec_controlledPhaseShift(qureg, idQubit1, idQubit2, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledPhaseShift(qureg, idQubit1+shift, idQubit2+shift, -angle);
    }
    
    qasm_recordControlledParamGate(qureg, GATE_PHASE_SHIFT, idQubit1, idQubit2, angle);
}

void multiControlledPhaseShift(QubitRegister qureg, int *controlQubits, int numControlQubits, REAL angle) {
    validateMultiControls(qureg, controlQubits, numControlQubits, __func__);
    
    statevec_multiControlledPhaseShift(qureg, controlQubits, numControlQubits, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(controlQubits, numControlQubits, shift);
        statevec_multiControlledPhaseShift(qureg, controlQubits, numControlQubits, angle);
        shiftIndices(controlQubits, numControlQubits, -shift);
    }
    
    qasm_recordMultiControlledParamGate(qureg, GATE_PHASE_SHIFT, controlQubits, numControlQubits-1, controlQubits[numControlQubits-1], angle);
}

void controlledNot(QubitRegister qureg, const int controlQubit, const int targetQubit) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledNot(qureg, controlQubit, targetQubit);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledNot(qureg, controlQubit+shift, targetQubit+shift);
    }
    
    qasm_recordControlledGate(qureg, GATE_SIGMA_X, controlQubit, targetQubit);
}

void controlledSigmaY(QubitRegister qureg, const int controlQubit, const int targetQubit) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledSigmaY(qureg, controlQubit, targetQubit);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledSigmaYConj(qureg, controlQubit+shift, targetQubit+shift);
    }
    
    qasm_recordControlledGate(qureg, GATE_SIGMA_Y, controlQubit, targetQubit);
}

void controlledPhaseFlip(QubitRegister qureg, const int idQubit1, const int idQubit2) {
    validateControlTarget(qureg, idQubit1, idQubit2, __func__);
    
    statevec_controlledPhaseFlip(qureg, idQubit1, idQubit2);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledPhaseFlip(qureg, idQubit1+shift, idQubit2+shift);
    }
    
    qasm_recordControlledGate(qureg, GATE_SIGMA_Z, idQubit1, idQubit2);
}

void multiControlledPhaseFlip(QubitRegister qureg, int *controlQubits, int numControlQubits) {
    validateMultiControls(qureg, controlQubits, numControlQubits, __func__);
    
    statevec_multiControlledPhaseFlip(qureg, controlQubits, numControlQubits);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(controlQubits, numControlQubits, shift);
        statevec_multiControlledPhaseFlip(qureg, controlQubits, numControlQubits);
        shiftIndices(controlQubits, numControlQubits, -shift);
    }
    
    qasm_recordMultiControlledGate(qureg, GATE_SIGMA_Z, controlQubits, numControlQubits-1, controlQubits[numControlQubits-1]);
}

void rotateAroundAxis(QubitRegister qureg, const int rotQubit, REAL angle, Vector axis) {
    validateTarget(qureg, rotQubit, __func__);
    validateVector(axis, __func__);
    
    statevec_rotateAroundAxis(qureg, rotQubit, angle, axis);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_rotateAroundAxisConj(qureg, rotQubit+shift, angle, axis);
    }
    
    qasm_recordAxisRotation(qureg, angle, axis, rotQubit);
}

void controlledRotateAroundAxis(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle, Vector axis) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    validateVector(axis, __func__);
    
    statevec_controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, axis);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledRotateAroundAxisConj(qureg, controlQubit+shift, targetQubit+shift, angle, axis);
    }
    
    qasm_recordControlledAxisRotation(qureg, angle, axis, controlQubit, targetQubit);
}

int getNumQubits(QubitRegister qureg) {
    return qureg.numQubitsRepresented;
}

int getNumAmps(QubitRegister qureg) {
    validateStateVecQureg(qureg, __func__);
    
    return qureg.numAmpsTotal;
}

REAL getRealAmpEl(QubitRegister qureg, long long int index) {
    validateStateVecQureg(qureg, __func__);
    validateStateIndex(qureg, index, __func__);
    
    return statevec_getRealAmpEl(qureg, index);
}

REAL getImagAmpEl(QubitRegister qureg, long long int index) {
    validateStateVecQureg(qureg, __func__);
    validateStateIndex(qureg, index, __func__);
    
    return statevec_getImagAmpEl(qureg, index);
}

REAL getProbEl(QubitRegister qureg, long long int index) {
    validateStateVecQureg(qureg, __func__);
    validateStateIndex(qureg, index, __func__);
    
    return statevec_getProbEl(qureg, index);
}

Complex getDensityAmplitude(QubitRegister qureg, long long int row, long long int col) {
    validateDensityMatrQureg(qureg, __func__);
    validateStateIndex(qureg, row, __func__);
    validateStateIndex(qureg, col, __func__);
    
    long long ind = row + col*(1LL << qureg.numQubitsRepresented);
    Complex amp;
    amp.real = statevec_getRealAmpEl(qureg, ind);
    amp.imag = statevec_getImagAmpEl(qureg, ind);
    return amp;
}

REAL calcTotalProbability(QubitRegister qureg) {
    if (qureg.isDensityMatrix)  
            return densmatr_calcTotalProbability(qureg);
        else
            return statevec_calcTotalProbability(qureg);
}

Complex calcInnerProduct(QubitRegister bra, QubitRegister ket) {
    validateStateVecQureg(bra, __func__);
    validateStateVecQureg(ket, __func__);
    validateMatchingQuregDims(bra, ket,  __func__);
    
    return statevec_calcInnerProduct(bra, ket);
}

REAL findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome) {
    validateTarget(qureg, measureQubit, __func__);
    validateOutcome(outcome, __func__);
    
    if (qureg.isDensityMatrix)
        return densmatr_findProbabilityOfOutcome(qureg, measureQubit, outcome);
    else
        return statevec_findProbabilityOfOutcome(qureg, measureQubit, outcome);
}

void cloneQubitRegister(QubitRegister targetQureg, QubitRegister copyQureg) {
    validateMatchingQuregTypes(targetQureg, copyQureg, __func__);
    validateMatchingQuregDims(targetQureg, copyQureg, __func__);
    
    statevec_cloneQubitRegister(targetQureg, copyQureg);
}

REAL calcPurity(QubitRegister qureg) {
    validateDensityMatrQureg(qureg, __func__);
    
    return densmatr_calcPurity(qureg);
}

REAL calcFidelity(QubitRegister qureg, QubitRegister pureState) {
    validateSecondQuregStateVec(pureState, __func__);
    validateMatchingQuregDims(qureg, pureState, __func__);
    
    if (qureg.isDensityMatrix)
        return densmatr_calcFidelity(qureg, pureState);
    else
        return statevec_calcFidelity(qureg, pureState);
}

void initPureState(QubitRegister qureg, QubitRegister pure) {
    validateSecondQuregStateVec(pure, __func__);
    validateMatchingQuregDims(qureg, pure, __func__);

    if (qureg.isDensityMatrix)
        densmatr_initPureState(qureg, pure);
    else
        statevec_cloneQubitRegister(qureg, pure);
    
    qasm_recordComment(qureg, "The register was initialised by a undisclosed given pure state.");
}

REAL collapseToOutcome(QubitRegister qureg, const int measureQubit, int outcome) {
    validateTarget(qureg, measureQubit, __func__);
    validateOutcome(outcome, __func__);
    
    REAL outcomeProb;
    if (qureg.isDensityMatrix) {
        outcomeProb = densmatr_findProbabilityOfOutcome(qureg, measureQubit, outcome);
        validateMeasurementProb(outcomeProb, __func__);
        densmatr_collapseToKnownProbOutcome(qureg, measureQubit, outcome, outcomeProb);
    } else {
        outcomeProb = statevec_findProbabilityOfOutcome(qureg, measureQubit, outcome);
        validateMeasurementProb(outcomeProb, __func__);
        statevec_collapseToKnownProbOutcome(qureg, measureQubit, outcome, outcomeProb);
    }
    
    qasm_recordMeasurement(qureg, measureQubit);
    return outcomeProb;
}

int measureWithStats(QubitRegister qureg, int measureQubit, REAL *outcomeProb) {
    validateTarget(qureg, measureQubit, __func__);

    int outcome;
    if (qureg.isDensityMatrix)
        outcome = densmatr_measureWithStats(qureg, measureQubit, outcomeProb);
    else
        outcome = statevec_measureWithStats(qureg, measureQubit, outcomeProb);
    
    qasm_recordMeasurement(qureg, measureQubit);
    return outcome;
}

int measure(QubitRegister qureg, int measureQubit) {
    validateTarget(qureg, measureQubit, __func__);
    
    int outcome;
    REAL discardedProb;
    if (qureg.isDensityMatrix)
        outcome = densmatr_measureWithStats(qureg, measureQubit, &discardedProb);
    else
        outcome = statevec_measureWithStats(qureg, measureQubit, &discardedProb);
    
    qasm_recordMeasurement(qureg, measureQubit);
    return outcome;
}

void addDensityMatrix(QubitRegister combineQureg, REAL otherProb, QubitRegister otherQureg) {
    validateDensityMatrQureg(combineQureg, __func__);
    validateDensityMatrQureg(otherQureg, __func__);
    validateMatchingQuregDims(combineQureg, otherQureg, __func__);
    validateProb(otherProb, __func__);
    
    densmatr_addDensityMatrix(combineQureg, otherProb, otherQureg);
}



/* new experimental dephasing functions */

// @TODO add to CPU local and distributed
void oneQubitDephase(QubitRegister qureg, const int targetQubit, REAL dephase) {
    validateDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, targetQubit, __func__);
    validateNoise(dephase, __func__);
    
    densmatr_oneQubitDephase(qureg, targetQubit, dephase);
}

// @TODO add to CPU local and distributed
void twoQubitDephase(QubitRegister qureg, const int qubit1, const int qubit2, REAL dephase) {
    validateDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, qubit1, __func__);
    validateTarget(qureg, qubit2, __func__);
    validateNoise(dephase, __func__);

    densmatr_twoQubitDephase(qureg, qubit1, qubit2, dephase);
}

// @TODO add to CPU local and distributed
void oneQubitDepolarise(QubitRegister qureg, const int targetQubit, REAL depolLevel) {
    validateDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, targetQubit, __func__);
    validateNoise(depolLevel, __func__);
    
    densmatr_oneQubitDepolarise(qureg, targetQubit, depolLevel);
}

// @TODO add to CPU local and distributed
void twoQubitDepolarise(QubitRegister qureg, const int qubit1, const int qubit2, REAL depolLevel) {
    validateDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, qubit1, __func__);
    validateTarget(qureg, qubit2, __func__);
    validateNoise(depolLevel, __func__);
    
    densmatr_twoQubitDepolarise(qureg, qubit1, qubit2, depolLevel);
}








/* debug and unit-testing functions */

int compareStates(QubitRegister qureg1, QubitRegister qureg2, REAL precision) {
    return statevec_compareStates(qureg1, qureg2, precision);
}

void initStateDebug(QubitRegister qureg) {
    statevec_initStateDebug(qureg);
}

void initStateFromSingleFile(QubitRegister *qureg, char filename[200], QuESTEnv env) {
    validateStateVecQureg(*qureg, __func__);

    int success = statevec_initStateFromSingleFile(qureg, filename, env);
    validateFileOpened(success, __func__);
}

void initStateOfSingleQubit(QubitRegister *qureg, int qubitId, int outcome) {
    return statevec_initStateOfSingleQubit(qureg, qubitId, outcome);
}

void reportStateToScreen(QubitRegister qureg, QuESTEnv env, int reportRank)  {
    statevec_reportStateToScreen(qureg, env, reportRank);
}


#ifdef __cplusplus
}
#endif