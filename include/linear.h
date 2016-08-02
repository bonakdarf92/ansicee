//
// Created by Farid Bonakdar on 30.06.16.
//

#ifndef ANSICEE_LINEAR_H
#define ANSICEE_LINEAR_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "horizontalModel.h"

#define K 1.0
#define R 3.0
#define S 12.0
#define I_MAX_A 10
#define I_MAX_B 1000
#define T_DELTA 0.1603
#define T_N_CONSTUP 0.0739
#define T_N_CONSTDN 0.6211
#define A0 10



/*
 * Kontrollieren ob alle Methoden Implementiert wurden
 * und pruefen ob evtl. etwas weiter modularisiert werden kann.
 */

// Test funktion


void initMatrix();

void matrixPresetting();

void getInputParameter();

void calculate_KS();

void calculate_KI();

void changing_engineSpeed(int n_updn);

void calculate_Cop();

void calculate_Dop();

double Matrix_KP(double bEnd, double Matrix_KS);

double EigenwertI(double EigenwertAI, int i, int a0);

void drehzahlAenderung(double Matrix[12][12]);


#endif //ANSICEE_LINEAR_H
