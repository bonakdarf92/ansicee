#include <stdio.h>
//#include "modelAntrieb.h"
//#include "DGL_Berechnen.h"
//#include "Subsysteme.h"
//#include "calcMatrix_A_B.h"
#include "horizontalModel.h"
#include "linear.h"
//#include "Lenkmotoren.h"
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
//#include "matrizeCalculator.h"
#include "InitTest.h"




int main() {

    double maxTime = 0;
    double minTime = 0.02;
    initializeVector();
    open_files();
    start_initializing();
    start_reading();
    initTest();
    initMatrix();
    matrixPresetting();
    size_t zaehler = 0;
    double gesamtzeit = 0;
    while (zaehler < 6000) {
        // Taking time for performance
        clock_t begin = clock();

        // This method combines all functions written in the File horizontalModel.c
        calculate_C_and_D(zaehler);

        getInputParameter();
        changing_engineSpeed(1);
        calculate_Cop();
        calculate_Dop();
        calculate_KS();
        calculate_KI(get_Matrix(4),scalar(1));
        calculate_Ai();
        calculate_EWI();
        matrix_Calculator_EWI();
        tune_matrix_EWI();
        calculate_KP(get_Matrix(4), scalar(2));
        calculate_AG();
        //calculate_EWG();
/*

        //Zum testen fuer Vectoren
        printf(" ug --> [");
        for (size_t l = 0; l < 9; l++) {
            printf("%g ", gsl_vector_get(getVector(2), l));
        }
        printf("]\n");
        printf(" xg --> [");
        for (size_t l = 0; l < 3; l++) {
            printf("%g ", gsl_vector_get(getVector(1), l));
        }
        printf("]\n");


        printf(" ug_alt --> [");
        for (size_t l = 0; l < 9; l++) {
            printf("%g ", gsl_vector_get(getVector(4), l));
        }
        printf("]\n");

        printf(" xg_alt --> [");
        for (size_t l = 0; l < 3; l++) {
            printf("%g ", gsl_vector_get(getVector(3), l));
        }
        printf("]\n");
        printf(" D --> [");

        for (size_t m = 0; m < 18; m++) {
            printf("%g ", gsl_vector_get(getMatrix(2), m));
        }
        printf("]\n");
        //saving_current_state();

*/

        //printf("%f\n%f\n", ax, ay);

        // Ende Zeitstop und Ausgabe der Zeit
        clock_t end = clock();
        double time = (double) 1000 * (end - begin) / CLOCKS_PER_SEC;
        if (maxTime <= time){
            maxTime = time;
        }
        //if (minTime>=time){
          //  minTime = time;
        //}
        //gesamtzeit += time;
        //printf("Berechnung dauert %.4f ms\n", time);
        zaehler++;
    }

    printf("Laengste Berechnung dauert %.4f ms\n", maxTime);
    //printf("Kuerzeste Berechnung dauert %.4f ms\n", minTime);

    return 0;
}
