//#include <stdio.h>
#include <stdint.h>
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
#include <sys/times.h>
#include <sys/time.h>
//#include "matrizeCalculator.h"
#include "InitTest.h"





int main() {

    double maxTime = 0;
    double minTime = 100;
    initializeVector();
    open_files();
    start_initializing();
    start_reading();
    initTest();
    initMatrix();
    matrixPresetting();
    size_t zaehler = 0;
    double gesamtzeit = 0;
    struct timeval tv;
    while (zaehler < 6000) {
        // Taking time for performance
        clock_t begin = clock();

        // This method combines all functions written in the File horizontalModel.c
        calculate_C_and_D(zaehler);

        // This method combines all function written in linear.c
        //calculating_PI_Controller();
        /*
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
        calculate_EWG();
        matrix_Calculator_EWG();
        //tune_matrix_EWG();
        tune_KP();
        */

/*
        //Zum testen fuer Vectoren
        printf(" delta_u --> [");
        for (size_t l = 0; l < 9; l++) {
            printf("%g ", gsl_vector_get(getVector(7), l));
        }
        printf("]\n");
        printf(" delta_x --> [");
        for (size_t l = 0; l < 3; l++) {
            printf("%g ", gsl_vector_get(getVector(6), l));
        }
        printf("]\n");
*/

        //printf("%f\n%f\n", ax, ay);

        // Ende Zeitstop und Ausgabe der Zeit
        clock_t end = clock();

        float time = (float) (end - begin);
        if (maxTime <= time){
            maxTime = time;
        }
        if (minTime>=time){
            minTime = time;
        }
        gesamtzeit += time;
        printf("Berechnung dauert %g s\n", time / CLOCKS_PER_SEC);
        zaehler++;
    }
    //struct tms tms1;
    //printf("Gesamtanzahl der Zyklen betraegt %f \n", (double )times(&tms1));
    printf("Laengste Berechnung dauert %g us\n", 1000000 * (maxTime/CLOCKS_PER_SEC));
    printf("Kuerzeste Berechnung dauert %g us\n", 1000000*(minTime/CLOCKS_PER_SEC));

    return 0;
}
