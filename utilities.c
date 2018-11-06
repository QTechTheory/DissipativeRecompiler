
#include <QuEST.h>
#include <stdio.h>
#include <stdlib.h>

#define DEFAULT_ARRAY_VALUE 0



/* multi-dimensional arrays */

double* createArray(int len) {
    
    double* arr = malloc(len * sizeof *arr);
    for (int i=0; i < len; i++)
        arr[i] = DEFAULT_ARRAY_VALUE;
    
    return arr;
}

double** createNestedArray(int outerLen, int innerLen) {
    
    double** arr = malloc(outerLen * sizeof *arr);
    for (int i=0; i < outerLen; i++) {
        arr[i] = malloc(innerLen * sizeof **arr);
        for (int j=0; j < innerLen; j++)
            arr[i][j] = DEFAULT_ARRAY_VALUE;
    }
    
    return arr;
}

void freeNestedArray(double** arr, int outerLen) {
    
    for (int i=0; i < outerLen; i++)
        free(arr[i]);
    free(arr);
}

void cloneArray(double* dest, double* source, int len) {
    
    for (int i=0; i < len; i++)
        dest[i] = source[i];
}



/* caching parameters */

double* loadParams(char* filename, int numParams) {
    
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf(
            "ERROR: loading params failed since file %s was not found! "    
            "Exiting...\n", filename);
        exit(1);
    }
    
    // skip descriptive header
    char line[2000];
    fgets(line, 2000*sizeof(char), file);
    
    // read each param on a newline
    double* params = malloc(numParams * sizeof *params);
    for (int p=0; p < numParams; p++) {
        
        // make sure we haven't tried to read too many
        if (fscanf(file, "%lf\n", &(params[p])) == EOF) {
            printf(
                "ERROR! Expected %d params, but file (%s) ended after reading %d\n",
                numParams, filename, p);
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }
    
    // make sure we didn't fewer than the file contained
    if (!feof(file)) {        
        printf(
            "ERROR! Expected %d params, but file (%s) contained more unread params!\n",
            numParams, filename);
        fclose(file);
        exit(EXIT_FAILURE);
    }
    
    // return the params which must be freed by caller
    fclose(file);
    return params;
}

void writeParams(char* filename, double* params, int numParams, char* header) {
    
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf(
            "ERROR: writing params failed since file %s could not be written to! "    
            "Exiting...\n", filename);
        exit(1);
    }
    
    
    fprintf(file, "%s\n", header);
    for (int p=0; p < numParams; p++)
        fprintf(file, "%.17g\n", params[p]);
    fclose(file);
}



/* caching wavefunctions */

void loadWavefunction(char* filename, QubitRegister qureg) {
    
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf(
            "ERROR: loading wavefunction failed since file %s was not found! "    
            "Exiting...\n", filename);
        exit(1);
    }
    
    // skip descriptive header
    char line[200];
    fgets(line, 200*sizeof(char), file);
    
    // read each param on a newline
    for (long long int i=0; i < qureg.numAmpsTotal; i++) {
        
        // first read each line to two buffers
        char strRe[50];
        char strIm[50];
                
        // make sure we haven't tried to read too many
        if (fscanf(file, "%s %s\n", strRe, strIm) == EOF) {
            printf(
                "ERROR! Expected %lld amplitudes, but file (%s) ended after reading %lld\n",
                qureg.numAmpsTotal, filename, i);
            fclose(file);
            exit(EXIT_FAILURE);
        }
        
        // update statevector with parsed buffers
        qureg.stateVec.real[i] = strtod(strRe, NULL);
        qureg.stateVec.imag[i] = strtod(strIm, NULL);
    }
    
    // make sure we didn't fewer than the file contained
    if (!feof(file)) {
        printf(
            "ERROR! Expected %lld amplitudes, but file (%s) contained more unread params!\n",
            qureg.numAmpsTotal, filename);
        fclose(file);
        exit(EXIT_FAILURE);
    }
    
    fclose(file);
}
    
void writeWavefunction(char* filename, QubitRegister qureg, char* header) {
    
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf(
            "ERROR: writing wavefunction failed: file %s could not be written to! "    
            "Exiting...\n", filename);
        exit(1);
    }
    
    
    fprintf(file, "%s\n", header);
    for (long long int i=0; i < qureg.numAmpsTotal; i++) {
        double ampRe = getRealAmpEl(qureg, i);
        double ampIm = getImagAmpEl(qureg, i);
        fprintf(file, "%.17g %.17g\n", ampRe, ampIm);
    }
    fclose(file);
}
