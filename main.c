#include <stdio.h>
//#include "modelAntrieb.h"
#include "DGL_Berechnen.h"
//#include "Subsysteme.h"
//#include "calcMatrix_A_B.h"
#include "horizontalModel.h"
#include "linear.h"
//#include "Lenkmotoren.h"
#include <time.h>
//#include "matrizeCalculator.h"
#include "InitTest.h"



int main() {

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
        clock_t begin = clock();
        testVector(zaehler);
        slipage();
        adma_velocity();
        slip();
        friction();
        double ax = Bewegungsgleichung_ax();
        double ay = Bewegungsgleichung_ay();
        gsl_vector *aufstand = AufstandsKraefte();
        RadKraefte();
        GierbewegungBerechnen();
        SystemmatrixBerechnen();
        deltasBerechnen();
        getInputParameter();
        calculate_Cop();
        //calculate_Dop();
        //calculate_KS();

        // Zum testen fuer Vectoren
        /*
        printf("[");
        for (int l = 0; l < 3; l++) {
            printf("%.3f ", gsl_vector_get(getVector(8), l));
        }
        printf("]");

        //Zum testen fuer Vectoren
        printf("\n[");
        for (int l = 0; l < 3; l++) {
            printf("%.3f ", gsl_vector_get(getVector(9), l));
        }
        printf("]\n");

        printf("[");
        for (int l = 0; l < 3; l++) {
            printf("%.3f ", gsl_vector_get(getVector(10), l));
        }
        printf("]\n");
        printf("[");

        for (int m = 0; m < 17; m++) {
            printf("%.2f ", gsl_matrix_get(getMatrix(1), m, 0));
        }
        printf("]\n");


        printf("[");

        for (int m = 0; m < 17; m++) {
            printf("%.2f ", gsl_matrix_get(getMatrix(2), m, 0));
        }
        printf("]\n");

        */

        printf("%f\n%f\n", ax, ay);
        // Ende Zeitstop und Ausgabe der Zeit
        clock_t end = clock();
        double time = (double) 1000 * (end - begin) / CLOCKS_PER_SEC;
        gesamtzeit += time;
        printf("Berechnung dauert %.4f ms\n", time);
        zaehler++;
    }
    return 0;
}
