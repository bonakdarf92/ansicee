#include <stdio.h>
//#include "modelAntrieb.h"
//#include "DGL_Berechnen.h"
//#include "Subsysteme.h"
//#include "calcMatrix_A_B.h"
#include "horizontalModel.h"
//#include "linear.h"
//#include "Lenkmotoren.h"
#include <time.h>
//#include "matrizeCalculator.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>



int main() {

    clock_t begin = clock();
    double daten[] = {4,-14,-2, -14, 65, 3, -2, 3, 3};
    gsl_matrix* mat = gsl_matrix_alloc(3,3);
    int k = 0;
    // Befuellen der Matrix m
    /*
    for (int i = 0; i <3 ; i++) {
        for (int j = 0; j < 3; j++) {
            gsl_matrix_set(mat, i, j, daten[k]);
            k++;
        }
    }
     */
    int i,j;
    /*
    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++) {
            printf("%.0f ", gsl_matrix_get(mat, i, j));
        }
        printf("\n");
    }
    */
    gsl_matrix_free(mat);
    initializeVector();
    testVector();
    schraeglaufwinkel();
    Adma_geschwindigkeit();
    schlupfBerechnung();
    ReibwertBerechnung();
    Bewegungsgleichung_ax();

    // Zum testen fuer Vectoren
    printf("[");
    for (int l = 0; l < 3; l++) {
        printf("%.3f ",gsl_vector_get(getVector(8),l));
    }
    printf("]");

    //Zum testen fuer Vectoren
    printf("\n[");
    for (int l = 0; l < 3; l++) {
        printf("%.3f ",gsl_vector_get(getVector(9),l));
    }
    printf("]\n");

    printf("[");
    for (int l = 0; l < 3; l++) {
        printf("%.3f ",gsl_vector_get(getVector(10),l));
    }
    printf("]\n");

    // Ende Zeitstop und Ausgabe der Zeit
    clock_t  end = clock();
    double time = (double) 1000*(end-begin)/CLOCKS_PER_SEC;
    printf("Berechnung dauert %.4f ms\n", time);
    return 0;
}
