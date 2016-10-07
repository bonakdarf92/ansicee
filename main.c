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
    clock_t alpha = clock();

    //if(RegelungOn==1)
    initializeVector();
    open_files(0);
    start_initializing(0);
    start_reading();
    initTest();
    initMatrix();
    matrixPresetting();
    size_t zaehler = 0;
    float timings[61000];
    while (zaehler < 61000) {
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

        zaehler++;
    }
    print_Timings(timings, sizeof(timings)/ sizeof(float));
    calculateCycleTime(timings, "MIN");
    calculateCycleTime(timings, "MAX");
    calculateCycleTime(timings, "Total");
    clock_t omega = clock();

    printf("Gesamte Berechnung dauert %g\n", (float) (omega - alpha) / CLOCKS_PER_SEC);

    //printf("Unterschied betraegt %f", calculate_difference(getVector(1),getVector(11)));
    return 0;
}
