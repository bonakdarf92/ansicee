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

    gsl_vector* correction;
    gsl_vector* tempSize3;
    gsl_vector* tempSize18;

    double korrekterAcc;
    tempSize3 = gsl_vector_alloc(3);
    tempSize18 = gsl_vector_alloc(18);
    correction = gsl_vector_alloc(3);
    initializeVector();
    open_files(0);
    start_initializing(0);
    start_reading();
    drehzahlTester();
    initTest();
    initMatrix();
    matrixPresetting();
    initCorrection();
    size_t zaehler = 0;
    float timings[61001];
    while (zaehler < 61001) {

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


        gsl_matrix_get_row(tempSize18, savingMatrix(2), zaehler);
        printf("Korrekter Wert: ");
        printer(NULL, tempSize18);
        correction = simple_difference(getMatrix(2), tempSize18);
        printf("Korrektur     : ");
        printer(NULL, correction);

        printf("Berechn. Wert : ");
        printer(NULL, getMatrix(2));

/*
        // Vergleich der Beschleunigungen
        korrekterAcc = gsl_vector_get(saving(15), zaehler);
        printf("Korrekter Wert: %f\n", korrekterAcc);
        printf("Berechner Wert: %f\n", gsl_vector_get(getVector(18), 2));
        printf("Differenz Wert: %f\n", gsl_vector_get(getVector(18), 2)-korrekterAcc);
*/

        //printf("                  n1       n2       n3       deldyn1   deldyn2   deldyn3   delta1   delta2    delta3\n");
        //printf(" Vektor ug    : ");
        //printer(NULL, getVector(2));

        printf("_________________________________\n");
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
