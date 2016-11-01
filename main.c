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
#include "calcMatrix_A_B.h"
#include "ErrorCorrection.h"
#define DEBUGGER 0




int main() {
    clock_t alpha = clock();

#if DEBUGGER == 1
    gsl_vector* correction;
    gsl_vector* tempSize3;
    gsl_vector* tempSize18;
    gsl_vector* tempSize7;
    gsl_vector* tempSize9;
    gsl_vector* tempSize4;

    double korrekterAcc;
    tempSize3 = gsl_vector_alloc(3);
    tempSize18 = gsl_vector_alloc(18);
    tempSize7 = gsl_vector_alloc(7);
    tempSize9 = gsl_vector_alloc(9);
    tempSize4 = gsl_vector_alloc(4);
    correction = gsl_vector_alloc(9);
#endif

    gsl_matrix* tempKI;
    tempKI = gsl_matrix_alloc(12, 3);
    gsl_matrix* tempProd1;
    tempProd1 = gsl_matrix_alloc(3, 3);
    initializeVector();
    open_files(0);
    start_initializing(0);
    start_reading();
    drehzahlTester();
    initTest();
    initMatrix();
    matrixPresetting();
    initCorrection();
    //double brain = 0;
    //size_t count;
    //double brain2 = 0;
    //double brain3 = 0;

    size_t zaehler = 0;
    float timings[61001];
    while (zaehler < 61001) {

        // Taking time for performance
        clock_t begin = clock();
        testVector(zaehler);
        calculate_C_and_D(zaehler);
        calculating_PI_Controller(zaehler);
        //create_data(zaehler, bigMatrices(2), tempKI);
        create_data(zaehler, bigMatrices(4), tempProd1);

        //printf("Korrekter KS ");
        //printer(tempKS, NULL);

        //printf("Prod berech. : ");
        //printer(get_Matrix(16), NULL);

        //printf("Prod korrekt : ");
        //printer(tempProd1, NULL);

        //printf("Difference : ");
        //printer(complex_difference(get_Matrix(16), tempProd1), NULL);


        //create_data(zaehler, bigMatrices(5), tempProd1);
        //printf("Inv.  Korr : ");
        //printer(tempProd1, NULL);
        //printf("Inv.  Ber. : ");
        //printer(get_Matrix(17), NULL);



#if DEBUGGER == 1
        gsl_matrix_get_row(tempSize18, savingMatrix(1), zaehler);
        printf("c Matrix Korr : ");
        printer(NULL, tempSize18);



        // Korrekte Zustandsvektor aus Simulink
        //gsl_matrix_get_row(tempSize3, savingMatrix(18), zaehler);
        //printf("Korr.  xg          : ");
        //printer(NULL, tempSize3);

        // Kopierte Werte aus C-Code
        //printf("Anf.    xg         : ");
        //printer(NULL, getVector(1));

        // Zustandsvektor aus vorherigen Iteration
        //printf("Anf. xg_alt        : ");
        //printer(NULL, getVector(3));

        // Differnz aus xg und xg_alt
        //printf("Anf. xg - xg_alt   : ");
        //printer(NULL, simple_difference(getVector(1), getVector(3)));

        // Korrekte Differenz aus Simulink
        printf("Korr. c Matrix : ");
        //gsl_matrix_get_row(tempSize3, savingMatrix(26), zaehler);
        printer(NULL, getMatrix(1));
        //printer(NULL,getVector(1));
#endif

        //printf(" \n +++++++ Vor Berechnung +++++++ \n \n");
        // This method combines all functions written in the File horizontalModel.c
        //calculate_C_and_D(zaehler);
        //printer(NULL, savingMatrix(1));

 /*

        printf("End    xg         : ");
        printer(NULL, getVector(1));

        printf("End    xg_alt     : ");
        printer(NULL, getVector(3));


        printf("End xg - xg_alt   : ");
        printer(NULL, simple_difference(getVector(1), getVector(3)));
        */

        //printf("Difference     : ");
        //printer(NULL, simple_difference(getMatrix(1), tempSize18));

        //printf("xg-xg_alt End : ");
        //printer(NULL, simple_difference(getVector(1), getVector(3)));

        // This method combines all function written in linear.c
        //calculating_PI_Controller();

        //printf("c Matrix  Ber.: ");
        //printer(get_Matrix(7), NULL);
        // Ende Zeitstop und Ausgabe der Zeit
        clock_t end = clock();

        // Speichern der Zeitpunkte in das Array
        timings[zaehler] = saveTiming(begin, end);



        //korrekterAcc = gsl_vector_get(saving(14), zaehler);
        //printf("%.10f\n", korrekterAcc);
        //printf("Berechnet     : ");
        //printer(NULL, correction);
        //printer(NULL, getVector(18));
        //correction = simple_difference(getVector(18), getVector(5));
        //printf("Berechn. Wert : %.11f\n", returnAcceleration(2));
        //printer(NULL, buildVector(1));
        //printf("acc_alt       : ");
        //printer(NULL, getVector(5));
        //printf("Endzustand xg : ");
        //printer(NULL, getVector(1));

        //printf("Endzu. xg_alt : ");
        //printer(NULL, getVector(3));

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

        //printf(" Vektor diff  : ");
        //printer(NULL, simple_difference(tempSize18, getMatrix(1)));
        //printer(get_Matrix(7), NULL);


        printf("\n ___________________________________________________\n \n %zu\n", zaehler+1);
        zaehler++;
    }
    //print_Timings(timings, sizeof(timings)/ sizeof(float));
    //calculateCycleTime(timings, "MIN");
    //calculateCycleTime(timings, "MAX");
    //calculateCycleTime(timings, "Total");
    clock_t omega = clock();
    printf("Gesamte Berechnung dauert %g\n", (float) (omega - alpha) / CLOCKS_PER_SEC);


    return 0;
}
