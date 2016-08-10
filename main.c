#include <stdio.h>
//#include "modelAntrieb.h"
#include "DGL_Berechnen.h"
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

    initializeVector();
    open_files();
    start_initializing();
    start_reading();
    initTest();
    initMatrix();
    matrixPresetting();
    size_t zaehler = 0;
    double gesamtzeit = 0;
    while (zaehler < 8) {
        printf("_______________ %zu  ______________\n",zaehler);
        clock_t begin = clock();
        testVector(zaehler);
        slipage();
        adma_velocity();
        slip();
        friction();
        double ax = Bewegungsgleichung_ax();
        double ay = Bewegungsgleichung_ay();
        AufstandsKraefte();
        RadKraefte();
        GierbewegungBerechnen();
        deltasBerechnen();
        SystemmatrixBerechnen();
        saving_current_state();
        //getInputParameter();
        //calculate_Cop();
        //calculate_Dop();
        //calculate_KS();
        //calculate_KI(get_Matrix(4),scalar(1));
        //calculate_Ai();


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
        saving_current_state();



        printf("%f\n%f\n", ax, ay);
        // Ende Zeitstop und Ausgabe der Zeit
        clock_t end = clock();
        double time = (double) 1000 * (end - begin) / CLOCKS_PER_SEC;
        gesamtzeit += time;
        printf("Berechnung dauert %.4f ms\n", time);
        zaehler++;
/*
        gsl_vector* ausgabe = get_Matrix(7);
        printf(" C_op -->\n");
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 12; j++) {
                printf("%f ", gsl_matrix_get(ausgabe, i,j));
            }
            printf("\n");
        }
        printf(" <--- C_op\n");
        */

        printf("++++++++++++++++++++++++++++++++++++\n");
    }


    return 0;
}
