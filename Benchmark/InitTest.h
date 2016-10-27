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


void open_files(size_t choice);

void start_initializing(size_t choice);

void start_reading(void);

gsl_vector* saving(size_t n);

gsl_matrix* savingMatrix(size_t n);

//gsl_matrix* bigMatrices(size_t n);

void create_data(size_t zaehler, gsl_matrix* matrix, gsl_matrix* matrixSaver);

void printer(gsl_matrix* matrix, gsl_vector* vector);

void calculateCycleTime (float timings[], char *auswahl);

float saveTiming(clock_t Anfang, clock_t Ende);

int compare_strings(char a[], char b[]);

void print_Timings(float timings[],size_t size);

void drehzahlTester(void);




#endif //ANSICEE_INITTEST_H
