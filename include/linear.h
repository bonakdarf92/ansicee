//
// Created by Farid Bonakdar on 30.06.16.
//

#ifndef ANSICEE_LINEAR_H
#define ANSICEE_LINEAR_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "horizontalModel.h"
// #include <gsl/gsl_blas.h>
#include <gsl/gsl_vector_complex.h>
#include <stdio.h>
#include <gsl/gsl_eigen.h>

#define K 1.0
#define R 3.0
#define S 12
#define I_MAX_A 10
#define I_MAX_B 1000
#define T_DELTA 0.1603
#define T_N_CONSTUP 0.0739
#define T_N_CONSTDN 0.6211




/*
 * Kontrollieren ob alle Methoden Implementiert wurden
 * und pruefen ob evtl. etwas weiter modularisiert werden kann.
 */

// Test funktion


void initMatrix();

void matrixPresetting();

gsl_matrix * get_Matrix(size_t n);

void getInputParameter();

void calculate_KS();

void calculate_KI(gsl_matrix* KS, double a0);

void calculate_Ai();

void calculate_EWI();

void changing_engineSpeed(int n_updn);

void calculate_Cop();

void calculate_Dop();

void matrix_Calculator_EWI();

void tune_matrix_EWI();

void calculate_KP(gsl_matrix* KS, double b0);

void calculating_AG();


#endif //ANSICEE_LINEAR_H
