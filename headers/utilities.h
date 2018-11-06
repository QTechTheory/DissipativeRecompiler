#ifndef UTILITIES_H_
#define UTILITIES_H_



/* multi-dimensional arrays */

double* createArray(int len);

double** createNestedArray(int outerLen, int innerLen);

void freeNestedArray(double** arr, int outerLen);

void cloneArray(double* dest, double* source, int len);


/* caching parameters */

double* loadParams(char* filename, int numParams);

void writeParams(char* filename, double* params, int numParams, char* header);



/* caching wavefunctions */

void loadWavefunction(char* filename, QubitRegister qureg);

void writeWavefunction(char* filename, QubitRegister qureg, char* header);



#endif // UTILITIES_H_