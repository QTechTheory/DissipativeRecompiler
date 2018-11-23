#ifndef CIRCUIT_LOADER_H_
#define CIRCUIT_LOADER_H_

#include <QuEST.h>

typedef enum {X, Y, Z, XX, YY, ZZ, Rx, Ry, Rz, RxRx, RyRy, RzRz, CTRL, cRx, cRy, cRz, Id} gateLabel;

typedef struct {
    
    int numQubits;
    int numGates;
    
    gateLabel* gateLabels;
    int* gateFirstQubits;
    int* gateSecondQubits;
    
} Circuit;

Circuit loadCircuit(char* filename);

void printGate(Circuit circuit, int g);

void printCircuit(Circuit circuit);

void applyCircuit(Circuit circuit, double* params, QubitRegister qureg);

void applyReversedCircuit(Circuit circuit, double* params, QubitRegister qureg);

void applyRandomisedCircuit(Circuit circuit, double* params, QubitRegister qureg);

void applyInverseCircuit(Circuit circuit, double* params, QubitRegister qureg);

void applyOrderedCircuit(Circuit circuit, double* params, int* ordering, QubitRegister qureg);

void applyReversedOrderedCircuit(Circuit circuit, double* params, int* ordering, QubitRegister qureg);

void applyGate(Circuit circuit, int gateNum, QubitRegister qureg);

void freeCircuit(Circuit circuit);

Circuit loadCircuitFromHamilFile(char* filename, int skipCoeff);

#endif // CIRCUIT_LOADER_H_