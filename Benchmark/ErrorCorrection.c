//
// Created by Farid Bonakdar on 05.10.16.
//

#include "ErrorCorrection.h"


gsl_rstat_workspace* workspace2;     // Workspace for statistical operations
gsl_vector* CMatrixTest;
gsl_vector* DMatrixTest;
gsl_vector* CMatrixRef;
gsl_vector* DMatrixRef;
gsl_vector* difference3;
gsl_vector* difference18;
double rmsErrorC [61001];
double rmsErrorD [61001];

/*
 * This function initializes the workspace for
 * all further calculations
 */
void initCorrection(){
    workspace2 = gsl_rstat_alloc();
    CMatrixRef = gsl_vector_alloc(18);
    DMatrixRef = gsl_vector_alloc(18);
    CMatrixTest = gsl_vector_alloc(18);
    DMatrixTest = gsl_vector_alloc(18);
    difference3 = gsl_vector_alloc(3);
    difference18 = gsl_vector_alloc(18);
}

/*
 * This method gets two vectors 1 and 2 and calculates the
 * root mean square of the difference of 1 and 2.
 * The lower the value the bigger the precision
 */
double calculate_difference(gsl_vector* one, gsl_vector* two){
    //initCorrection();

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

/*
 * This Method gets two input vectors, one calculated from the horizontal model
 * and a reference one
 */
gsl_vector* simple_difference(gsl_vector* calculated, gsl_vector* reference){
    gsl_vector_memcpy(difference3, calculated);
    gsl_vector_sub(difference3, reference);
    return difference3;
}

/*
 * This Method picks the Vectors C and D from the C-Code
 * and the reference Vectors from the Simulink Model
 */
void storeCurrentInformation(size_t cyc){
    gsl_matrix_get_row(CMatrixRef, savingMatrix(1), cyc);
    gsl_matrix_get_row(DMatrixRef, savingMatrix(2), cyc);
    CMatrixTest = getMatrix(1);
    DMatrixTest = getMatrix(2);
}


void calculateCorrection(size_t cyc){
    //initCorrection();
    storeCurrentInformation(cyc);
    //rmsErrorC[cyc] = calculate_difference(getMatrix(1), CMatrixRef);
    //rmsErrorD[cyc] = calculate_difference(getMatrix(2), DMatrixRef);
    difference3 = simple_difference(getMatrix(1), CMatrixRef);
    size_t j;
    for (j = 0; j < 18 ; ++j) {
        printf("%.4f ",gsl_vector_get(difference3, j));
    }
    printf("\n");
}

gsl_vector* returnReference(size_t n){
    switch (n){
        case 1:
            return CMatrixRef;
        case 2:
            return DMatrixRef;
    }
}

