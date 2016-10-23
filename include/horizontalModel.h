//
// Created by Farid Bonakdar on 30.06.16.
//

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdlib.h>

#ifndef ANSICEE_HORIZONTALMODEL_H
#define ANSICEE_HORIZONTALMODEL_H

#define KMU 5.0                     // linearized mu-slip about 5
#define C_a 9095.00000              // slippage
#define RADIUS 0.15                 // radius of tire
#define MASSE 1200.0000             // total mass kg originally 308.52
#define HCG 0.48303                 // height of the center of gravity m  -- CarMaker self calculated // 0.55m platform
#define LAENGE 2.28                 // equilateral triangle setup m initially 2.54
#define GRAVY 9.81                  // gravity constant
#define THETA 944.8465962           // Angle Theta
#define SQRT3 1.7320508075          // square root of 3

// If you want to Debug and see whats happening change the value to 1
#define DEBUGGER 0

#if DEBUGGER == 1
    double ax1;
    double ax2;
    double ax3;
    double ax4;
    double ax5;
    double ax6;
    double ax7;
    double ay1;
    double ay2;
    double ay3;
    double ay4;
#endif

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


void initializeVector(void);

gsl_vector * getVector(size_t n);

gsl_vector * getMatrix(size_t n);

void testVector(size_t n);

void initTest(void);

void slipage(void);

void adma_velocity(void);

void slip(void);

void friction(void);

void acceleration_ax(void);

void acceleration_ay(void);

void contactForce(void);

void wheelForce(void);

void yawrateCalculator(void);

void systemMatrixCalculator(void);

void deltasBerechnen(void);

void saving_current_state(void);

void calculate_C_and_D(size_t cyc);

double returnAcceleration(size_t n);

gsl_vector* buildVector(size_t choice);

#endif //ANSICEE_HORIZONTALMODEL_H
