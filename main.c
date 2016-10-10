//#include <stdio.h>
#include <stdint.h>
#include "horizontalModel.h"
#include "linear.h"
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <sys/times.h>
#include <sys/time.h>
#include "InitTest.h"
#include "ErrorCorrection.h"




int main() {

    gsl_vector* tempnachher;
    tempnachher = gsl_vector_alloc(3);


    gsl_vector* tempjetzt;
    tempjetzt = gsl_vector_alloc(3);
    clock_t alpha = clock();

    //if(RegelungOn==1)
    initializeVector();
    open_files(0);
    start_initializing(0);
    start_reading();
    initTest();
    initMatrix();
    matrixPresetting();
    initCorrection();
    size_t zaehler = 0;
    float timings[1001];
    while (zaehler < 1001) {

        //tempnachher = getVector(3);
        //gsl_vector_sub(temp, getVector(3));
        //printf("xg_alt.: %zu", zaehler);
        //printer(NULL, tempnachher);
        saving_current_state();
        storeCurrentInformation(zaehler);


        // Taking time for performance
        clock_t begin = clock();

        // This method combines all functions written in the File horizontalModel.c
        calculate_C_and_D(zaehler);

        // This method combines all function written in linear.c
        calculating_PI_Controller();
        //Zum testen fuer Vectoren
        //printer(get_Matrix(7), NULL);

        // Ende Zeitstop und Ausgabe der Zeit
        clock_t end = clock();

        // Speichern der Zeitpunkte in das Array
        timings[zaehler] = saveTiming(begin, end);

        //gsl_matrix_get_row(temp, savingMatrix(1), zaehler);
        //calculateCorrection(zaehler);

        //gsl_vector_sub(temp, getVector(3));
        //printf("xg_mom.: %zu",zaehler);
        //printer(NULL, getMatrix(1));

        //tempjetzt = getVector(1);
        //printf("xg_alt.: %zu",zaehler);
        //printer(NULL, getVector(3));
        //printer(NULL, tempnachher);



        //tempjetzt = getVector(1);
        printf("Cref.: %zu", zaehler);
        printer(NULL, returnReference(1));
        printf("Cakt.: %zu", zaehler);
        printer(NULL,getMatrix(1));
        zaehler++;
    }
    //print_Timings(timings, sizeof(timings)/ sizeof(float));
    //calculateCycleTime(timings, "MIN");
    //calculateCycleTime(timings, "MAX");
    //calculateCycleTime(timings, "Total");
    clock_t omega = clock();

    printf("Gesamte Berechnung dauert %g\n", (float) (omega - alpha) / CLOCKS_PER_SEC);

    //printf("Unterschied betraegt %f", calculate_difference(getVector(1),getVector(11)));
    return 0;
}
