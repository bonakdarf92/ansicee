//
// Created by Farid Bonakdar on 30.06.16.
//

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#ifndef ANSICEE_HORIZONTALMODEL_H
#define ANSICEE_HORIZONTALMODEL_H

//TODO Anpassen des Infotextes
/*
 * This Model is transcribed from Jan's Horizontal Model.
 * All Methods and Calculation formula are taken from his
 * Matlab-Code or his Masterthesis.
 * The main purpose of this c-file is to generate the System
 * Matrices C and D.
 * Some thought should be spent on the model of capsulising
 * the Information. A
 * For further Information take a look in the readme, Bachelor-
 * thesis of Farid Bonakdar or Masterthesis of Jan Steier
 */

void initializeVector();

gsl_vector * getVector(int n);

gsl_matrix * getMatrix(int n);

void testVector();

void slipage();

void adma_velocity();

void slip();

void friction();

double Bewegungsgleichung_ax();

double Bewegungsgleichung_ay();

gsl_vector * AufstandsKraefte();

void RadKraefte();

void GierbewegungBerechnen();

void SystemmatrixBerechnen();

void deltasBerechnen();

#endif //ANSICEE_HORIZONTALMODEL_H
