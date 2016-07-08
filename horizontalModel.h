//
// Created by Farid Bonakdar on 30.06.16.
//

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#ifndef ANSICEE_HORIZONTALMODEL_H
#define ANSICEE_HORIZONTALMODEL_H

//TODO Anpassen des Infotextes
/*
 * This Model is transcripted from Jan's Horizontalmodel.
 * All Methods and Calculation formula are taken from his
 * Matlab-Code or his Masterthesis.
 * For futher Information take a look in the readme, Bachelor-
 * thesis of Farid Bonakdar or Masterthesis of Jan Steier
 */

void initializeVector();

gsl_vector * getVector(int n);

gsl_matrix * getMatrix(int n);

void testVector();

void schraeglaufwinkel();

void Adma_geschwindigkeit();

void schlupfBerechnung();

void ReibwertBerechnung();

double Bewegungsgleichung_ax();

double Bewegungsgleichung_ay();

gsl_vector *  AufstandsKraefte();

void RadKraefte();

void GierbewegungBerechnen();

void SystemmatrixBerechnen();

void deltasBerechnen();

#endif //ANSICEE_HORIZONTALMODEL_H
