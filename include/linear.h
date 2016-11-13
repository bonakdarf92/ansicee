//
// Created by Farid Bonakdar on 30.06.16.
//

#ifndef ANSICEE_LINEAR_H
#define ANSICEE_LINEAR_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "horizontalModel.h"
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_linalg.h>
#include <Lapacke/lapacke.h>
#include <Lapacke/lapacke_utils.h>
#include <Lapacke/lapacke_config.h>


#define KONSTANTE 1.0
#define I_MAX_A 10
#define I_MAX_B 1000
#define T_DELTA 0.160256410256410
#define T_N_CONSTUP 0.073909830007391
#define T_N_CONSTDN 0.621118012422360




/*
 * Kontrollieren ob alle Methoden Implementiert wurden
 * und pruefen ob evtl. etwas weiter modularisiert werden kann.
 */


gsl_matrix* ASystem;        // Declaration of Matrix A
gsl_matrix* BSystem;        // Declaration of Matrix B
gsl_matrix* A_Inv;          // Declaration of Matrix A inverse
gsl_matrix* KI;             // Declaration of Matrix KI
gsl_matrix* KS;             // Declaration of Matrix KS
gsl_matrix* KP;
gsl_matrix* Ai;             // Declaration of Matrix A iterator
gsl_matrix* Ag;
gsl_matrix_complex* EW_I;           // Declaration of Matrix Eigenvalue for Integrator
gsl_matrix* EW_G;           // Declaration of Matrix Eigenvalue for G??
gsl_matrix* EW_I1;          // Declaration of Matrix Eigenvalue for inner calculations
gsl_matrix* EW_G1;
gsl_matrix* EWI_output;     // Declaration of output Matrix EWI
gsl_matrix* EWG_output;

gsl_vector* C_in;           // Declaration of input Matrix C
gsl_vector* D_in;           // Declaration of input Matrix D

gsl_matrix* C_op;           // Declaration of operation Matrix C
gsl_matrix* D_op;           // Declaration of operation Matrix D
gsl_matrix* temp1;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temp2;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temp3;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temp4;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temp5;
gsl_matrix* temp6;
gsl_matrix* temporary;
gsl_matrix* Inverse;
gsl_matrix* eye;            // Declaration of identity matrix for calculations KP
double a_min = 0.001;       // Declaration and initialization of tuning factor a_min
double A0 = 10.0;
double b_min = 0.00001;     // Declaration and initialization of tuning factor b_min
double B0 = 10.0;             // TODO check which value B0 has

gsl_eigen_nonsymm_workspace* workspace;     // Declaration of workspace for eigenvalue calculation
gsl_vector_complex* eigenvalue;             // Declaration of vector for eigenvalues
gsl_vector_complex* eigenvalue2;            // Declaration of vector for eigenvalues
gsl_vector_complex* eigenvalue3;            // Declaration of vector for eigenvalues
gsl_vector_complex* eigenvalue4;            // Declaration of vector for eigenvalues
gsl_matrix* dmd;
gsl_matrix* threetimes3;
gsl_matrix* linalgtest;
gsl_matrix* KI2;
gsl_matrix* KI3;
gsl_matrix* KI4;

gsl_vector_view a;

void initMatrix(void);

void matrixPresetting(void);

gsl_matrix * get_Matrix(size_t n);

void getInputParameter(void);

void calculate_KS(void);

void calculate_KI(gsl_matrix* KS, double a0);

void calculate_Ai(void);

void calculate_EWI(void);

void changing_engineSpeed(gsl_matrix* n_updn, size_t zaehler);

void calculate_Cop(void);

void calculate_Dop(void);

void matrix_Calculator_EWI(void);

void tune_matrix_EWI(void);

void calculate_KP(gsl_matrix* KS, double b0);

void calculate_AG(void);

void calculate_EWG(void);

void matrix_Calculator_EWG(void);

void tune_matrix_EWG(void);

void tune_KP(void);

double scalar(size_t n);

void calculating_PI_Controller(size_t zaehler);

gsl_vector* returnCinDin(size_t n);

void update_BSystem(size_t zaehler);

void rechanging_engineSpeed(void);


#endif //ANSICEE_LINEAR_H
