//
// Created by Farid Bonakdar on 04.08.16.
//

#ifndef ANSICEE_INITTEST_H
#define ANSICEE_INITTEST_H

//#include <stdio.h>
#include <wchar.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


void open_files(void);

void start_initializing(void);

void start_reading(void);

gsl_vector* saving(size_t n);

void create_data(void);

void printer(gsl_matrix* matrix, gsl_vector* vector);

void calculateCycleTime (float timings[], char *auswahl);

float saveTiming(clock_t Anfang, clock_t Ende);

int compare_strings(char a[], char b[]);




#endif //ANSICEE_INITTEST_H
