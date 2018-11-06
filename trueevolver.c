
#include "trueevolver.h"
#include "hamiltonianloader.h"

#include <QuEST.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

/* signatures */

void populatePropogator(TrueEvolEnv evEnv, Hamiltonian hamil, double timeStep);

/* definitions */

TrueEvolEnv initTrueEvolEnv(Hamiltonian hamil, double timeStep, QuESTEnv qEnv) {
    
    TrueEvolEnv evEnv;
    evEnv.numQubits = hamil.numQubits;
    evEnv.timeStep = timeStep;
    evEnv.dim = 2LL << (evEnv.numQubits - 1);
    evEnv.copyStateInVec = gsl_vector_complex_alloc(evEnv.dim);
    evEnv.copyStateOutVec = gsl_vector_complex_alloc(evEnv.dim);  
    evEnv.propogator = gsl_matrix_complex_alloc(evEnv.dim, evEnv.dim);
    populatePropogator(evEnv, hamil, timeStep);
    
    return evEnv;
}

void closeTrueEvolEnv(TrueEvolEnv evEnv, QuESTEnv qEnv) {
    gsl_vector_complex_free(evEnv.copyStateInVec);
    gsl_vector_complex_free(evEnv.copyStateOutVec);
    gsl_matrix_complex_free(evEnv.propogator);
}

/** Taken from:
 * https://stackoverflow.com/questions/10049160/complex-matrix-exponential-in-c
 */
void my_gsl_complex_matrix_exponential(
    gsl_matrix_complex *eA, gsl_matrix_complex *A, int dimx)
{
    int j,k=0;
    gsl_complex temp;
    gsl_matrix *matreal =gsl_matrix_alloc(2*dimx,2*dimx);
    gsl_matrix *expmatreal =gsl_matrix_alloc(2*dimx,2*dimx);
    
    //Converting the complex matrix into real one using A=[Areal, Aimag;-Aimag,Areal]
    for (j = 0; j < dimx;j++)
        for (k = 0; k < dimx;k++)
        {
            temp=gsl_matrix_complex_get(A,j,k);
            gsl_matrix_set(matreal,j,k,GSL_REAL(temp));
            gsl_matrix_set(matreal,dimx+j,dimx+k,GSL_REAL(temp));
            gsl_matrix_set(matreal,j,dimx+k,GSL_IMAG(temp));
            gsl_matrix_set(matreal,dimx+j,k,-GSL_IMAG(temp));
        }

    // prec mode: https://github.com/ampl/gsl/blob/master/gsl_mode.h
    // uses method 3 of https://www.cs.cornell.edu/cv/ResearchPDF/19ways+.pdf
    gsl_linalg_exponential_ss(matreal,expmatreal,GSL_PREC_DOUBLE);

    double realp;
    double imagp;
    for (j = 0; j < dimx;j++)
        for (k = 0; k < dimx;k++)
        {
            realp=gsl_matrix_get(expmatreal,j,k);
            imagp=gsl_matrix_get(expmatreal,j,dimx+k);
            gsl_matrix_complex_set(eA,j,k,gsl_complex_rect(realp,imagp));
        }
    gsl_matrix_free(matreal);
    gsl_matrix_free(expmatreal);
}

void populatePropogator(TrueEvolEnv evEnv, Hamiltonian hamil, double timeStep) {
    
    // create -iHt matrix
    gsl_matrix_complex* powerMatr = gsl_matrix_complex_alloc(evEnv.dim, evEnv.dim);
    gsl_matrix_complex_memcpy(powerMatr, hamil.matrix);
    gsl_matrix_complex_scale(powerMatr, gsl_complex_rect(0,-timeStep));
    
    // create propogator = exp(-iHt)
    my_gsl_complex_matrix_exponential(evEnv.propogator, powerMatr, evEnv.dim);
    
    //@debug: printing
    //gsl_matrix_complex_fprintf(stdout, evEnv.propogator, "%g");
    
    // free -iHT matrix
    gsl_matrix_complex_free(powerMatr);
}

void evolveTrueState(QubitRegister qureg, Hamiltonian hamil, TrueEvolEnv evEnv) {
    
    // copying to a vector is unnecessarily expensive, but eh
    for (int i=0; i < evEnv.dim; i++)
        gsl_vector_complex_set(evEnv.copyStateInVec, i, 
            gsl_complex_rect(getRealAmpEl(qureg, i), getImagAmpEl(qureg, i)));
    
    // apply propogator (outVec= 1*propogator*outVec + 0*outVec)
    gsl_blas_zgemv(
        CblasNoTrans, gsl_complex_rect(1,0), evEnv.propogator, evEnv.copyStateInVec, 
        gsl_complex_rect(0,0), evEnv.copyStateOutVec);
        
    // copy copyStateOutVec into qureg
    for (int i=0; i < evEnv.dim; i++) {
        gsl_complex amp = gsl_vector_complex_get(evEnv.copyStateOutVec, i);
        qureg.stateVec.real[i] = GSL_REAL(amp);
        qureg.stateVec.imag[i] = GSL_IMAG(amp);
    }
}

