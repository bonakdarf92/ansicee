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
    correction = gsl_vector_alloc(9);
    initializeVector();
    open_files(0);
    start_initializing(0);
    start_reading();
    drehzahlTester();
    initTest();
    initMatrix();
    matrixPresetting();
    initCorrection();
    double brain = 0;
    size_t count;
    double brain2 = 0;
    double brain3 = 0;

    size_t zaehler = 0;
    float timings[61000];
    while (zaehler < 61000) {

        storeCurrentInformation(zaehler);

        //saving_current_state();
        // Taking time for performance
        clock_t begin = clock();
        //testVector(zaehler);
        //printf("vor der Berech: ");
        //printer(NULL, getVector(19));

        // This method combines all functions written in the File horizontalModel.c
        calculate_C_and_D(zaehler);

        // This method combines all function written in linear.c
        //calculating_PI_Controller();

        // Ende Zeitstop und Ausgabe der Zeit
        clock_t end = clock();

        // Speichern der Zeitpunkte in das Array
        timings[zaehler] = saveTiming(begin, end);


        gsl_matrix_get_row(tempSize3, savingMatrix(9), zaehler);
        printf("Korrekter Wert: ");
        printer(NULL, tempSize3);
        //korrekterAcc = gsl_vector_get(saving(14), zaehler);
        //printf("%.10f\n", korrekterAcc);
        //correction = simple_difference(getVector(20), tempSize3);
        printf("Berechnet     : ");
        //printer(NULL, correction);

        //printf("Berechn. Wert : %.10f\n", returnAcceleration(1));
        printer(NULL, getVector(19));
        //printf("%f\n", gsl_vector_get(getVector(2), 0));
        //printer(NULL, getVector(19));
/*
        // Vergleich der Beschleunigungen
        korrekterAcc = gsl_vector_get(saving(15), zaehler);
        printf("Korrekter Wert: %f\n", korrekterAcc);
        printf("Berechner Wert: %f\n", gsl_vector_get(getVector(18), 2));
*/      //printf("Differenz Wert: %.10f\n", fabs(returnAcceleration(1)-korrekterAcc));

        //if (fabs(returnAcceleration(1)-korrekterAcc) > brain){
        //    brain = fabs(returnAcceleration(1)-korrekterAcc);
        //    count = zaehler;
        //    brain2 = korrekterAcc;
        //    brain3 = returnAcceleration(1);
        //}

        //printf("                  n1       n2       n3       deldyn1   deldyn2   deldyn3   delta1   delta2    delta3\n");
        //printf(" Vektor diff  : ");
        //printer(NULL, simple_difference(tempSize18, getMatrix(1)));
        //printer(NULL, getVector(1));


        printf("_________________________________\n");
        zaehler++;
    }
    //print_Timings(timings, sizeof(timings)/ sizeof(float));
    //calculateCycleTime(timings, "MIN");
    //calculateCycleTime(timings, "MAX");
    //calculateCycleTime(timings, "Total");
    clock_t omega = clock();
    printf("Gesamte Berechnung dauert %g\n", (float) (omega - alpha) / CLOCKS_PER_SEC);
    //printf("Groesste Abweichung %f\nStelle %d\nKorrekter Wert %f \nBerechneter Wert %f\n", brain, count, brain2, brain3);

    //printf("Unterschied betraegt %f", calculate_difference(getVector(1),getVector(11)));
    return 0;
}
