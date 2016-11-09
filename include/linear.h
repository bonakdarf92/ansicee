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

// Test funktion


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
