//
// Created by Farid Bonakdar on 05.10.16.
//

#include "ErrorCorrection.h"


gsl_rstat_workspace* workspace2;     // Workspace for statistical operations
gsl_vector* cmatrix;
gsl_vector* dmatrix;


/*
 * This function initializes the workspace for
 * all further calculations
 */
void initCorrection(){
    workspace2 = gsl_rstat_alloc();
    //cmatrix = gsl_matrix_alloc()
}

/*
 * This method gets two vectors 1 and 2 and calculates the
 * root mean square of the difference of 1 and 2.
 * The lower the value the bigger the precision
 */
double calculate_difference(gsl_vector* one, gsl_vector* two){
    initCorrection();

    // Initialize output value, rang of both vectors and its sizes
    double rms;
    size_t a, b, i;
    a = one->size;
    b = two->size;

    // Declaration and initialization of current vector for subtraction
    gsl_vector* temp;
    temp = gsl_vector_alloc(a);

    // If vectors have different sizes print this error message and return negative value
    if(a != b || a == 0 || b == 0){
        rms = -1.0;
        printf("Fehler!!! Vektoren besiten nicht gleiche Laenge oder Null");
    }

    else {
        // Copy values of vector one into temp and subtract
        gsl_vector_memcpy(temp, one);
        gsl_vector_sub(temp, two);

        // Add each vector value to workspace for statistic calculation
        for (i = 0; i < a ; i++) {
            gsl_rstat_add(gsl_vector_get(temp,i),workspace2);
        }

        // Return root mean square of difference vector
        rms = gsl_rstat_rms(workspace2);
    }

    return rms;
}
