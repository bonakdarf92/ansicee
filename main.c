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

    if(RegelungOn==1)
    initializeVector();
    open_files();
    start_initializing();
    start_reading();
    initTest();
    initMatrix();
    matrixPresetting();
    size_t zaehler = 0;
    float timinings[6000];
    while (zaehler < 6000) {
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
        timinings[zaehler] = saveTiming(begin, end);

        zaehler++;
    }
    size_t j;
    for (j = 0;j<6000;j++){
        printf(" Zeit : %g\n",timinings[j]);
    }
    calculateCycleTime(timinings, "MIN");
    calculateCycleTime(timinings, "MAX");
    calculateCycleTime(timinings, "Total");
    return 0;
}
