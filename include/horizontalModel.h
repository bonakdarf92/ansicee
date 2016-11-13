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

/*
 * Declaration of System matrices and inner System vector for further use in methods.
 * Stored here to get access to each scalar, vector and matrix.
 */
gsl_vector* xg;             // Vector xg
gsl_vector* ug;             // Vector ug
gsl_vector* xg_alt;         // Vector xg_alt
gsl_vector* ug_alt;         // Vector ug_Alt
gsl_vector* acc_alt;        // Vector acc_alt
gsl_vector* delta_x;        // Vector delta_x
gsl_vector* delta_u;        // Vector delta_u
gsl_vector* alpha_r;        // Vector alpha_r
gsl_vector* alpha_x;        // Vector alpha_x
gsl_vector* alpha_y;        // Vector alpha_y
gsl_vector* C;              // Vector C
gsl_vector* D;              // Vector D
double v_1;                 // Velocity No.1
double v_2;                 // Velocity No.2
double v_3;                 // Velocity No.3
double beta_1;              // Angle beta_1
double beta_2;              // Angle beta_2
double beta_3;              // Angle beta_3
gsl_vector* beta;           // Vector beta
gsl_vector* v;              // Vector v
gsl_vector* vr;             // Vector vr
gsl_vector* sr;             // Vector sr
gsl_vector* mu;             // Vector mu
gsl_vector* mux;            // Vector mu
gsl_vector* muy;            // Vector mu
gsl_vector* Fz;             // Force Fz
gsl_vector* Fx;             // Force Fx
gsl_vector* Fy;             // Force Fy
gsl_vector* acc;            // Vector acc
gsl_matrix* test_ug;        // Matrix for testing and simulation ug
gsl_matrix* test_xg;        // Matrix for testing and simulation xg
double ax;
double ay;
double psi_pp;
gsl_vector* difference_xg_alt;
gsl_vector* difference_ug_alt;
gsl_vector* xg_alt_loop;
gsl_vector* ug_alt_loop;

/*
 * This Model is translated from Jan's Horizontal Model.
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
