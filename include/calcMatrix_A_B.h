//
// Created by Farid Bonakdar on 30.06.16.
//

#ifndef ANSICEE_CALCMATRIX_A_B_H
#define ANSICEE_CALCMATRIX_A_B_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>


void invertMatrixA(gsl_matrix* matrix, gsl_matrix* output);

void invertMatrixB(double Matrix[12][12]);

double eigenWerte();

double determinante(gsl_matrix* matrix);

void matrixMulti(double a[12][12], double b [12][12]);

void matrixAdd(double a[12][12], double b[12][12]);

void matrixSub(double a[12][12], double b[12][12]);

#endif //ANSICEE_CALCMATRIX_A_B_H
