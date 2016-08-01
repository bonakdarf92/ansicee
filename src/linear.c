//
// Created by Farid Bonakdar on 30.06.16.
//
#include "linear.h"


gsl_matrix* A;              // Declaration of Matrix A
gsl_matrix* B;              // Declaration of Matrix B
gsl_matrix* A_Inv;          // Declaration of Matrix A inverse
gsl_matrix* KI;             // Declaration of Matrix KI
gsl_matrix* KS;             // Declaration of Matrix KS
gsl_matrix* Ai;             // Declaration of Matrix A iterator
gsl_matrix* EW_I;           // Declaration of Matrix Eigenvalue for Integrator
gsl_matrix* EW_G;           // Declaration of Matrix Eigenvalue for G??
gsl_matrix* C_in;           // Declaration of input Matrix C
gsl_matrix* D_in;           // Declaration of input Matrix D



void initMatrix(){
    A = gsl_matrix_alloc(12,12);
    B = gsl_matrix_alloc(12,12);
    A_Inv = gsl_matrix_alloc(12,12);
    KS = gsl_matrix_alloc(12,3);
    // TODO initialize all matrices


    gsl_matrix_set_zero(A);
    gsl_matrix_set_zero(B);
    gsl_matrix_set_zero(A_Inv);

}

void matrixPresetting(){

    /*
     * This function puts all manually calculated indexes into the zero matrices A, B and A inverse
     * to get a better Overview look at the shape of the all matrices
     * Every matrix has a size of 12 x 12
     * index numbers begin with 0 not 1 in the code !!!!
     *
     *                                      Shape of these Matrices
     *              +---------------------------------------------------------------------------+
     *              | 1_1    1_2   1_3   1_4   1_5   1_6   1_7   1_8   1_9   1_10   1_11   1_12 |
     *              | 2_1    2_2   2_3   2_4   2_5   2_6   2_7   2_8   2_9   2_10   2_11   2_12 |
     *              | 3_1    3_2   3_3   3_4   3_5   3_6   3_7   3_8   3_9   3_10   3_11   3_12 |
     *              | 4_1    4_2   4_3   4_4   4_5   4_6   4_7   4_8   4_9   4_10   4_11   4_12 |
     *              | 5_1    5_2   5_3   5_4   5_5   5_6   5_7   5_8   5_9   5_10   5_11   5_12 |
     *              | 6_1    6_2   6_3   6_4   6_5   6_6   6_7   6_8   6_9   6_10   6_11   6_12 |
     *    Matrix =  | 7_1    7_2   7_3   7_4   7_5   7_6   7_7   7_8   7_9   7_10   7_11   7_12 |
     *              | 8_1    8_2   8_3   8_4   8_5   8_6   8_7   8_8   8_9   8_10   8_11   8_12 |
     *              | 9_1    9_2   9_3   9_4   9_5   9_6   9_7   9_8   9_9   9_10   9_11   9_12 |
     *              | 10_1  10_2  10_3  10_4  10_5  10_6  10_7  10_8  10_9  10_10  10_11  10_12 |
     *              | 11_1  11_2  11_3  11_4  11_5  11_6  11_7  11_8  11_9  11_10  11_11  11_12 |
     *              | 12_1  12_2  12_3  12_4  12_5  12_6  12_7  12_8  12_9  12_10  12_11  12_12 |
     *              +---------------------------------------------------------------------------+
     */


    // Setting all the indexes for Matrix A
    gsl_matrix_set(A, 0, 0, -7.2384);           // Setting value in Matrix A_1_1
    gsl_matrix_set(A, 0, 1, -38.9376);          // Setting value in Matrix A_1_2
    gsl_matrix_set(A, 1, 0, 1);                 // Setting value in Matrix A_2_1
    gsl_matrix_set(A, 2, 2, -7.2384);           // Setting value in Matrix A_3_3
    gsl_matrix_set(A, 2, 3, -38.9376);          // Setting value in Matrix A_3_4
    gsl_matrix_set(A, 3, 2, 1);                 // Setting value in Matrix A_4_3
    gsl_matrix_set(A, 4, 4, -7.2384);           // Setting value in Matrix A_5_5
    gsl_matrix_set(A, 4, 5, -38.9376);          // Setting value in Matrix A_5_6
    gsl_matrix_set(A, 5, 4, 1);                 // Setting value in Matrix A_6_5
    gsl_matrix_set(A, 6, 6, -7.0356);           // Setting value in Matrix A_7_7
    gsl_matrix_set(A, 6, 7, -183.0609);         // Setting value in Matrix A_7_8
    gsl_matrix_set(A, 7, 6, 1);                 // Setting value in Matrix A_8_7
    gsl_matrix_set(A, 8, 8, -7.0356);           // Setting value in Matrix A_9_9
    gsl_matrix_set(A, 8, 9, -183.0609);         // Setting value in Matrix A_9_10
    gsl_matrix_set(A, 9, 8, 1);                 // Setting value in Matrix A_10_9
    gsl_matrix_set(A, 10, 10, -7.0356);         // Setting value in Matrix A_11_11
    gsl_matrix_set(A, 10, 11, -183.0609);       // Setting value in Matrix A_11_12
    gsl_matrix_set(A, 11, 10, 1);               // Setting value in Matrix A_12_11


    // Setting all the indexes for Matrix B
    gsl_matrix_set(B, 0, 0, K / (T_DELTA*T_DELTA) );
    gsl_matrix_set(B, 2, 1, K / (T_DELTA*T_DELTA) );
    gsl_matrix_set(B, 4, 3, K / (T_DELTA*T_DELTA));
    gsl_matrix_set(B, 0, 0, -7.2384);
    gsl_matrix_set(B, 0,0,-7.2384);
    gsl_matrix_set(B, 0,0,-7.2384);
    gsl_matrix_set(B, 0,0,-7.2384);
    gsl_matrix_set(B, 0,0,-7.2384);
    gsl_matrix_set(B, 0,0,-7.2384);
    gsl_matrix_set(B, 0,0,-7.2384);
    gsl_matrix_set(B, 0,0,-7.2384);


    // Setting all the indexes for Matrix A inverse
    gsl_matrix_set(A_Inv, 0, 1, 1);             // Setting value in Matrix A_inv 1_1
    gsl_matrix_set(A_Inv, 1, 0, -0.0257);       // Setting value in Matrix A_inv 2_1
    gsl_matrix_set(A_Inv, 1, 1, 0.1859);        // Setting value in Matrix A_inv 2_2
    gsl_matrix_set(A_Inv, 2, 3, 1);             // Setting value in Matrix A_inv 3_4
    gsl_matrix_set(A_Inv, 3, 2, -0.0257);       // Setting value in Matrix A_inv 4_3
    gsl_matrix_set(A_Inv, 3, 3, 0.1859);        // Setting value in Matrix A_inv 4_4
    gsl_matrix_set(A_Inv, 4, 5, 1);             // Setting value in Matrix A_inv 5_6
    gsl_matrix_set(A_Inv, 5, 4, -0.0257);       // Setting value in Matrix A_inv 6_5
    gsl_matrix_set(A_Inv, 5, 5, 0.1859);        // Setting value in Matrix A_inv 6_6
    gsl_matrix_set(A_Inv, 6, 7, 1);             // Setting value in Matrix A_inv 7_8
    gsl_matrix_set(A_Inv, 7, 6, -0.0055);       // Setting value in Matrix A_inv 8_7
    gsl_matrix_set(A_Inv, 7, 7, -0.0384);       // Setting value in Matrix A_inv 8_8
    gsl_matrix_set(A_Inv, 8, 9, 1);             // Setting value in Matrix A_inv 9_10
    gsl_matrix_set(A_Inv, 9, 8, -0.0055);       // Setting value in Matrix A_inv 10_9
    gsl_matrix_set(A_Inv, 9, 9, -0.0384);       // Setting value in Matrix A_inv 10_10
    gsl_matrix_set(A_Inv, 10, 11, 1);           // Setting value in Matrix A_inv 11_12
    gsl_matrix_set(A_Inv, 11, 10, -0.0055);     // Setting value in Matrix A_inv 12_11
    gsl_matrix_set(A_Inv, 11, 11, -0.0384);     // Setting value in Matrix A_inv 12_12

}

void getInputParameter(){
    C_in = getMatrix(1);
    D_in = getMatrix(2);

}

//TODO Methode zur Berechnung der K_p Matrix
/*
 * Eingabeparameter unter anderem Tuningfaktor bEnd
 * und die Matrix K_s.
 * Auf weitere Input parameter pruefen
 */
double Matrix_KP(double bEnd, double Matrix_KS) {
    return 0;
}

//TODO Methode zur Berechnung der Eigenewert des I-Reglers

/*
 * Vorerst die Laufvariablen aus der Matlab datei uebernommen.
 * Bei Gelegenheit ueberpruefen ob diese benoetigt werden und
 * ob vereinfachung vorgenommen werden koennen.
 */

double EigenwertI(double EigenwertAI, int i, int a0) {
    return 0;
}


//TODO Methode fuer Matrixmanipulation wenn Drehzahl sich veraendert
/*
 * Laut Code von Jan muessen die Systemmatrizen A und Ainverse bei
 * Richtungsaenderung oder grossen Drehzahlaenderungen veraendert werden.
 * Fall abfrage aus Matlab Modell entnehmen
 */
void drehzahlAenderung(double (*Matrix)[12]) {

}

