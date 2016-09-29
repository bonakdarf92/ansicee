//
// Created by Farid Bonakdar on 04.08.16.
//

#ifndef ANSICEE_INITTEST_H
#define ANSICEE_INITTEST_H

//#include <stdio.h>
#include <wchar.h>
#include <gsl/gsl_vector.h>


void open_files(void);

void start_initializing(void);

void start_reading(void);

gsl_vector* saving(size_t n);

void create_data(void);

void printer(void);




#endif //ANSICEE_INITTEST_H
