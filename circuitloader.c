#include "circuitloader.h"

#include <QuEST.h>
#include <stdio.h>
#include <stdlib.h>

gateLabel getLabel(char c) {
    
    if (c == 'x')
        return Rx;
    if (c == 'y')
        return Ry;
    if (c == 'z')
        return Rz;
    if (c == 'X')
        return X;
    if (c == 'Y')
        return Y;
    if (c == 'Z')
        return Z;
    if (c == 'c')
        return CTRL;
    if (c == 'I')
        return Id;
    
    printf("ERROR! Unrecognised gate in circuit: %c\nExiting...\n", c);
    exit(1);
}

gateLabel getDoubleLabel(gateLabel l) {
    if (l == Rx)
        return RxRx;
    if (l == Ry)
        return RyRy;
    if (l == Rz)
        return RzRz;
    if (l == X)
        return XX;
    if (l == Y)
        return YY;
    if (l == Z)
        return ZZ;
    
    printf("INTERNAL ERROR! Unrecognised circuit gate enum! Exiting...\n");
    exit(1);
}

gateLabel getControlLabel(gateLabel l) {
    if (l == Rx)
        return cRx;
    if (l == Ry)
        return cRy;
    if (l == Rz)
        return cRz;
    
    printf("INTERNAL ERROR! Unrecognised circuit gate enum! Exiting...\n");
    exit(1);
}

char* getString(gateLabel l) {
    if (l == Rx)
        return "Rx";
    if (l == Ry)
        return "Ry";
    if (l == Rz)
        return "Rz";
    if (l == RxRx)
        return "RxRx";
    if (l == RyRy)
        return "RyRy";
    if (l == RzRz)
        return "RzRz";
    if (l == X)
        return "X";
    if (l == Y)
        return "Y";
    if (l == Z)
        return "Z";
    if (l == XX)
        return "XX";
    if (l == YY)
        return "YY";
    if (l == ZZ)
        return "ZZ";
    if (l == cRx)
        return "cRx";
    if (l == cRy)
        return "cRy";
    if (l == cRz)
        return "cRz";
    if (l == Id)
        return "I";

    printf("INTERNAL ERROR! Unrecognised circuit gate enum! Exiting...\n");
    exit(1);
}

Circuit loadCircuit(char *filename) {
    return loadCircuitFromHamilFile(filename, 0);
}

Circuit loadCircuitFromHamilFile(char* filename, int skipCoeff) {
    
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
		printf("ERROR: file (%s) not found!\n", filename);
        exit(1);
    }
    char ch;
    
    // count the number of gates
    int numGates = 1;
    while ((ch=getc(file)) != EOF)
        if (ch == '\n')
            numGates++;

    // allocate data structures
    gateLabel* gateLabels = malloc(numGates * sizeof *gateLabels);
    int* gateFirstQubits = malloc(numGates * sizeof *gateFirstQubits);
    int* gateSecondQubits = malloc(numGates * sizeof *gateSecondQubits);
    
    // find gates
    rewind(file);
    for (int g=0; g < numGates; g++) {
        
        // skip coeff (used in Hamiltonian spec files)
        if (skipCoeff) {
            double dummy;
            fscanf(file, "%lf ", &dummy);
        }
        
        // record single-qubit gate
        ch = getc(file);
        gateLabels[g] = getLabel(ch);
        fscanf(file, "%d", &(gateFirstQubits[g]));
                
        ch = getc(file);
        if (ch == '\n' || ch == EOF) {
            gateSecondQubits[g] = -1;
        } else {
            
            // besides controls, only identical ops are supported
            if (gateLabels[g] != CTRL && gateLabels[g] != getLabel(ch)) {
                printf(
                    "ERROR: Unrecognised circuit gate %s%c! Exiting...\n",
                    getString(gateLabels[g]), ch);
                exit(1);
            }
            
            // collet double-qubit gate info
            if (gateLabels[g] == CTRL)
                gateLabels[g] = getControlLabel(getLabel(ch));
            else
                gateLabels[g] = getDoubleLabel(gateLabels[g]);
            fscanf(file, "%d", &(gateSecondQubits[g]));
            getc(file); // burn \n
        }        
    }
    fclose(file);
    
    // shift qubit indices, find number of qubits
    int numQubits = 0;
    for (int g=0; g < numGates; g++) {
        
        if (gateFirstQubits[g] > numQubits)
            numQubits = gateFirstQubits[g];
        if (gateSecondQubits[g] > numQubits)
            numQubits = gateSecondQubits[g];
        
        gateFirstQubits[g] -= 1;
        if (gateSecondQubits[g] != -1)
            gateSecondQubits[g] -= 1;
    }
    
    // create circuit
    Circuit circuit;
    circuit.numQubits = numQubits;
    circuit.numGates = numGates;
    circuit.gateLabels = gateLabels;
    circuit.gateFirstQubits = gateFirstQubits;
    circuit.gateSecondQubits = gateSecondQubits;
    
    return circuit;
}

void printGate(Circuit circuit, int g) {
    
    gateLabel l = circuit.gateLabels[g];
    printf("%s on %d", getString(l), circuit.gateFirstQubits[g]);
    if (circuit.gateSecondQubits[g] != -1)
        printf(" and %d", circuit.gateSecondQubits[g]);
    printf("\n");
}

void printCircuit(Circuit circuit) {
    
    printf(
        "Circuit has %d gates acting on %d qubits. They are:\n", 
        circuit.numGates, circuit.numQubits
    );
    
    for (int g=0; g < circuit.numGates; g++)
        printGate(circuit, g);
}

void freeCircuit(Circuit circuit) {
    
    free(circuit.gateLabels);
    free(circuit.gateFirstQubits);
    free(circuit.gateSecondQubits);
}

/**
 * mode: 0=forward, 1=backward, 2=random, 3=inverse, 4=ordered, 5=backward ordered
 */
void applyCircuitInner(Circuit circuit, double* params, QubitRegister qureg, const int mode, int* ordering) {
    
    int gateWasApplied[circuit.numGates];
    for (int g=0; g < circuit.numGates; g++)
        gateWasApplied[g] = 0;
    
    // assumes unique parameter in every gate
    for (int t=0; t < circuit.numGates; t++) {

        // decide which gate to apply
        int g;
        if (mode == 0)
             g = t;
        else if (mode == 1 || mode == 3)
            g = circuit.numGates - 1 - t;
        else if (mode == 2) {
            g = (rand() * (long long int) circuit.numGates) / RAND_MAX;
            while (gateWasApplied[g]) g = (g + 1) % circuit.numGates;
            gateWasApplied[g] = 1;
        }
        else if (mode == 4)
            g = ordering[t];
        else if (mode == 5)
            g = ordering[circuit.numGates - 1 - t];
        else {
            printf(
                "INTERNAL ERROR: unrecognised mode %d passed to applyCircuitInner! Exiting...\n",
                mode);
            exit(1);
        }
        
        gateLabel gate = circuit.gateLabels[g];
        int qb1 = circuit.gateFirstQubits[g];
        int qb2 = circuit.gateSecondQubits[g];
        double param = params[g];
        
        // reverse the param to apply the inverse
        if (mode == 3)
            param *= -1;
        
        if (gate == Rx) {
            
            rotateX(qureg, qb1, param);   
        }
        else if (gate == Ry) {
            
            rotateY(qureg, qb1, param);
        }
        else if (gate == Rz) {
            
            rotateZ(qureg, qb1, param);
        }
        else if (gate == cRx) {
            
            controlledRotateX(qureg, qb1, qb2, param);
        }
        else if (gate == cRy) {
            
            controlledRotateY(qureg, qb1, qb2, param);
        }
        else if (gate == cRz) {
            
            controlledRotateZ(qureg, qb1, qb2, param);
        }
        else if (gate == RxRx) {
            
            controlledNot(qureg, qb1, qb2);
            rotateX(qureg, qb1, param);
            controlledNot(qureg, qb1, qb2);
        }
        else if (gate == RyRy) {
            
            controlledSigmaY(qureg, qb1, qb2);
            rotateY(qureg, qb1, param);
            controlledSigmaY(qureg, qb1, qb2);
        }
        else if (gate == RzRz) {
            
            controlledNot(qureg, qb1, qb2);
            rotateZ(qureg, qb2, param);
            controlledNot(qureg, qb1, qb2);
            
        } else {
            printf(
                "ERROR! Non-parameterised gate %s in ansatz circuit! Exiting...\n",
                getString(gate));
            exit(1);
        }
    }
}

void applyCircuit(Circuit circuit, double* params, QubitRegister qureg) {
    
    applyCircuitInner(circuit, params, qureg, 0, NULL);
}

void applyReversedCircuit(Circuit circuit, double* params, QubitRegister qureg) {
    
    applyCircuitInner(circuit, params, qureg, 1, NULL);
}

void applyRandomisedCircuit(Circuit circuit, double* params, QubitRegister qureg) {
    
    applyCircuitInner(circuit, params, qureg, 2, NULL);
}

void applyInverseCircuit(Circuit circuit, double* params, QubitRegister qureg) {
    
    applyCircuitInner(circuit, params, qureg, 3, NULL);
}

void applyOrderedCircuit(Circuit circuit, double* params, int* ordering, QubitRegister qureg) {
    
    applyCircuitInner(circuit, params, qureg, 4, ordering);
}

void applyReversedOrderedCircuit(Circuit circuit, double* params, int* ordering, QubitRegister qureg) {

    applyCircuitInner(circuit, params, qureg, 5, ordering);
}


/** only applies Pauli gates */
void applyGate(Circuit circuit, int gateNum, QubitRegister qureg) {
    
    gateLabel gate = circuit.gateLabels[gateNum];
    int qb1 = circuit.gateFirstQubits[gateNum];
    int qb2 = circuit.gateSecondQubits[gateNum];
    
    if (gate == X)
        sigmaX(qureg, qb1);
    else if (gate == Y)
        sigmaY(qureg, qb1);
    else if (gate == Z)
        sigmaZ(qureg, qb1);
    else if (gate == XX) {
        sigmaX(qureg, qb1);
        sigmaX(qureg, qb2);
    }
    else if (gate == YY) {
        sigmaY(qureg, qb1);
        sigmaY(qureg, qb2);
    }
    else if (gate == ZZ) {
        sigmaZ(qureg, qb1);
        sigmaZ(qureg, qb2);
    }
    else {
        printf(
            "ERROR! Non-Pauli gate %s in Hamiltonian term! Exiting...\n", 
            getString(gate));
        exit(1);
    }
}
