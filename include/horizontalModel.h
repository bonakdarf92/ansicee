//
// Created by Farid Bonakdar on 30.06.16.
//

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#ifndef ANSICEE_HORIZONTALMODEL_H
#define ANSICEE_HORIZONTALMODEL_H

#define KMU 5                         // linearized mu-slip about 5
#define C_a 9095                    // slippage
#define RADIUS 0.15                 // radius of tire
#define MASSE 1200                      // total mass kg originally 308.52
#define HCG 0.48303                 // height of the center of gravity m  -- CarMaker self calculated // 0.55m platform
#define LAENGE 2.28                      // equilateral triangle setup m initially 2.54
#define GRAVY 9.81                      // gravity constant
#define THETA 944.8465962           // Angle Theta
#define SQRT3 1.73205               // square root of 3

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

gsl_vector * getVector(size_t n);

gsl_vector * getMatrix(size_t n);

void testVector(size_t n);

void initTest();

void slipage();

void adma_velocity();

void slip();

void friction();

double Bewegungsgleichung_ax();

double Bewegungsgleichung_ay();

void AufstandsKraefte();

void RadKraefte();

void GierbewegungBerechnen();

void SystemmatrixBerechnen();

void deltasBerechnen();

void saving_current_state();

#endif //ANSICEE_HORIZONTALMODEL_H
