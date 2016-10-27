//
// Created by Farid Bonakdar on 05.10.16.
//

#ifndef ANSICEE_ERRORCORRECTION_H
#define ANSICEE_ERRORCORRECTION_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rstat.h>
#include "horizontalModel.h"
#include "InitTest.h"

void initCorrection();

double calculate_difference(gsl_vector* one, gsl_vector* two);

gsl_vector* simple_difference(gsl_vector* calculated, gsl_vector* refernce);

void storeCurrentInformation(size_t cyc);

void calculateCorrection(size_t cyc);

gsl_vector* returnReference(size_t n);

gsl_matrix* complex_difference(gsl_matrix* a, gsl_matrix* b);


#endif //ANSICEE_ERRORCORRECTION_H
