//
// Created by Farid Bonakdar on 04.08.16.
//

#ifndef ANSICEE_INITTEST_H
#define ANSICEE_INITTEST_H

#include <wchar.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "ErrorCorrection.h"
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex.h>


#if DEBUGGER == 1
FILE* Cmat;
FILE* Dmat;
FILE* alphar;
FILE* alphax;
FILE* alphay;
FILE* axFile;
FILE* ayFile;
FILE* betaFile;
FILE* FzFile;
FILE* FxFile;
FILE* FyFile;
FILE* muFile;
FILE* muxFile;
FILE* muyFile;
FILE* psiPPFile;
FILE* vFile;
FILE* vrFile;
FILE* stepUFile;
FILE* stepXFile;
#endif
FILE* DrehzahlenFile;
FILE* statLenkFile;
FILE* dynLenkFile;
FILE* xgFile;
FILE* ugFile;
FILE* accFile;
FILE* accAltFile;
#if DEBUGGER == 1
FILE* axtermeFile;
FILE* aytermeFile;
FILE* acc2File;
FILE* xg_diff_loopFile;
FILE* xg_alt_diffFile;
FILE* ug_diff_loopFile;
FILE* ug_alt_diffFile;
FILE* acc_diff_loopFile;
#endif
//FILE* EWgFile;
//FILE* EWiFile;
//FILE* EWgAnfangFile;
//FILE* EWgwhile1File;
//FILE* EWiAnfangFile;
//FILE* EWiwhile1File;
FILE* KiFile;
FILE* KiAnfangFile;
FILE* KiEndeFile;
FILE* KpFile;
FILE* KpAnfangFile;
FILE* KpEndeFile;
FILE* KptransFile;
FILE* KsOut1File;
FILE* AgAnfangFile;
FILE* AiOutFile;
FILE* n_updnFile;
FILE* T_nFile;
FILE* prod1File;
FILE* prod2File;

gsl_vector* n1;
gsl_vector* n2;
gsl_vector* n3;
gsl_vector* delta1;
gsl_vector* delta2;
gsl_vector* delta3;
gsl_vector* deltaDyn1;
gsl_vector* deltaDyn2;
gsl_vector* deltaDyn3;
gsl_vector* velocity_x;
gsl_vector* velocity_y;
gsl_vector* yawrate;
gsl_matrix* cMatrix;
gsl_matrix* dMatrix;
gsl_vector* a_x;
gsl_vector* a_y;
gsl_matrix* alphaRMatrix;
gsl_matrix* alphaXMatrix;
gsl_matrix* alphaYMatrix;
gsl_matrix* betaMatrix;
gsl_matrix* FxMatrix;
gsl_matrix* FyMatrix;
gsl_matrix* FzMatrix;
gsl_matrix* muMatrix;
gsl_matrix* muXMatrix;
gsl_matrix* muYMatrix;
gsl_vector* psiPPVec;
gsl_matrix* stepXMatrix;
gsl_matrix* stepUMatrix;
gsl_matrix* vMatrix;
gsl_matrix* vrMatrix;
gsl_matrix* DrehzahlenMat;
gsl_matrix* statLenkMatrix;
gsl_matrix* dynLenkMatrix;
gsl_matrix* xgMatrix;
gsl_matrix* ugMatrix;
gsl_matrix* accMatrix;
gsl_matrix* accAltMatrix;
gsl_matrix* axtermeMatrix;
gsl_matrix* aytermeMatrix;
gsl_matrix* acc2Matrix;
gsl_matrix* xg_diff_loopMatrix;
gsl_matrix* xg_alt_diffMatrix;
gsl_matrix* ug_diff_loopMatrix;
gsl_matrix* ug_alt_diffMatrix;
gsl_matrix* acc_diff_loopMatrix;
float gesamtzeit = 0.0;
//gsl_matrix_complex* EWgMatrix;
//gsl_matrix_complex* EWiMatrix;
//gsl_matrix_complex* EWgAnfangMatrix;
//gsl_matrix_complex* EWgwhile1Matrix;
//gsl_matrix_complex* EWiAnfangMatrix;
//gsl_matrix_complex* EWiwhile1Matrix;
gsl_matrix* KiMatrix;
gsl_matrix* KiAnfangMatrix;
gsl_matrix* KiEndeMatrix;
gsl_matrix* KpMatrix;
gsl_matrix* KpAnfangMatrix;
gsl_matrix* KpEndeMatrix;
gsl_matrix* KptransMatrix;
gsl_matrix* KsOut1Matrix;
gsl_matrix* AgAnfangMatrix;
gsl_matrix* AiOutMatrix;
gsl_vector* tempCol;
gsl_matrix* n_updnMatrix;
gsl_matrix* T_nMatrix;
gsl_matrix* prod1Matrix;
gsl_matrix* prod2Matrix;


void open_files(size_t choice);

void start_initializing(size_t choice);

void start_reading(void);

gsl_vector* saving(size_t n);

gsl_matrix* savingMatrix(size_t n);

gsl_matrix* bigMatrices(size_t n);

void create_data(size_t zaehler, gsl_matrix* matrix, gsl_matrix* matrixSaver);

void printer(gsl_matrix* matrix, gsl_vector* vector);

void calculateCycleTime (float timings[], char *auswahl);

float saveTiming(clock_t Anfang, clock_t Ende);

int compare_strings(char a[], char b[]);

void print_Timings(float timings[],size_t size);

void drehzahlTester(void);




#endif //ANSICEE_INITTEST_H
