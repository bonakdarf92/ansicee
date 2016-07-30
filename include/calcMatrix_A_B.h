//
// Created by Farid Bonakdar on 30.06.16.
//

#ifndef ANSICEE_CALCMATRIX_A_B_H
#define ANSICEE_CALCMATRIX_A_B_H


void invertMatrixA(double Matrix[12][12]);

void invertMatrixB(double Matrix[12][12]);

double eigenWerte();

double determinante();

void matrixMulti(double a[12][12], double b [12][12]);

void matrixAdd(double a[12][12], double b[12][12]);

void matrixSub(double a[12][12], double b[12][12]);

#endif //ANSICEE_CALCMATRIX_A_B_H
