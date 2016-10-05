//
// Created by Farid Bonakdar on 05.10.16.
//

#ifndef ANSICEE_ERRORCORRECTION_H
#define ANSICEE_ERRORCORRECTION_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rstat.h>

void initCorrection();

double calculate_difference(gsl_vector* one, gsl_vector* two);


#endif //ANSICEE_ERRORCORRECTION_H
