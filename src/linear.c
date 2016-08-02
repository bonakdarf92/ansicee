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

// TODO check what kind of data C_in and D_in are: Vector or Matrix
gsl_matrix* C_in;           // Declaration of input Matrix C
gsl_matrix* D_in;           // Declaration of input Matrix D

gsl_matrix* C_op;           // Declaration of operation Matrix C
gsl_matrix* D_op;           // Declaration of operation Matrix D
gsl_matrix* temp1;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temp2;          // Declaration of temporary Matrix for Calculations
double a_min = 0.001;       // Declaration and initialization of tuning factor a_min


/*
 * Initialing all matrices and setting up some matrix with zero for easy setup in further code
 * These steps are necessary because the calculations are made every iteration
 */
void initMatrix(){
    A = gsl_matrix_alloc(12,12);
    B = gsl_matrix_alloc(12,12);
    A_Inv = gsl_matrix_alloc(12,12);
    KS = gsl_matrix_alloc(12,3);
    C_op = gsl_matrix_alloc(12,3);
    D_op = gsl_matrix_alloc(12,3);
    temp1 = gsl_matrix_alloc(12,12);
    temp2 = gsl_matrix_alloc(12,3);
    // TODO initialize all matrices


    gsl_matrix_set_zero(A);
    gsl_matrix_set_zero(B);
    gsl_matrix_set_zero(A_Inv);
    gsl_matrix_set_zero(C_op);
    gsl_matrix_set_zero(D_op);

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

    //TODO T_n1, T_n2 und T_n3
    // Setting all the indexes for Matrix B
    gsl_matrix_set(B, 0, 0, K / (T_DELTA*T_DELTA) );
    gsl_matrix_set(B, 2, 1, K / (T_DELTA*T_DELTA) );
    gsl_matrix_set(B, 4, 2, K / (T_DELTA*T_DELTA));
    gsl_matrix_set(B, 6, 3, K / (T_N_CONSTUP*T_N_CONSTDN));     // TODO T_n_1
    gsl_matrix_set(B, 8, 4, K / (T_N_CONSTDN*T_N_CONSTDN));     // TODO T_n_2
    gsl_matrix_set(B, 10, 5, K / (T_N_CONSTDN*T_N_CONSTDN));    // TODO T_n_3


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


/*
 * This Method calls the function getMatrix from horizontal model and
 * store its returns in Matrix C_in and D_in.
 */

void getInputParameter(){
    C_in = getMatrix(1);
    D_in = getMatrix(2);
}


/*
 * This method calculates the matrix KS
 * The formula is taken from Jan's linearized model
 *
 *                         product 2
 *                      |-----+------|
 * --->         KS = D - C * A_inv * B
 * --->                     |----+---|
 * --->                      product 1
 * --->
 * --->             |--------+-------|
 *                      difference 3
 *
 * For performance purposes temporary matrices are generated to calculate
 * product 1, product 2 and at least difference 3
 */
void calculate_KS(){
    KS = D_op;
    temp1 = A_Inv;
    temp2 = C_op;
    gsl_matrix_mul_elements(temp1, B);
    gsl_matrix_alloc(temp2, temp1);
    gsl_matrix_sub(D_op, temp2);

}
// TODO might be a BOTTLENECK !!!
/*
 * This method calculates the matrix KI
 * The formula is taken from Jan's linearized model
 *
 *                         division 2
 *                    |------+-------|
 * --->     KS = A0 * (KS' / (KS*KS'))
 * --->                     |----+---|
 * --->                      product 1
 * --->
 * --->         |----------+---------|
 *                      const 3
 *
 * For performance purposes temporary matrices are generated to calculate
 * product 1, product 2 and at least difference 3
 */


void calculate_KI(){
    KI = D_op;

    temp1 = gsl_matrix_transpose(KS);
    temp2 = C_op;
    gsl_matrix_mul_elements(temp1, B);
    gsl_matrix_alloc(temp2, temp1);
    gsl_matrix_sub(D_op, temp2);

}
// TODO might be a BOTTLENECK !!!


/*
 * This method changes the matrix indexes for direction change in engine speed
 * For example the direction of movement converts from forward to reverse
 * Then the changes have to be made in System matrix A and A_inverse
 */

// TODO check what information is given in n_up down and how the differential observation takes place

void changing_engineSpeed(int n_updn){

    switch (n_updn) {
        case 1:
            gsl_matrix_set(A, 6, 6, -1.9964);
            gsl_matrix_set(A, 6, 7, -2.5921);
            gsl_matrix_set(A_Inv, 7, 6, -0.3858);
            gsl_matrix_set(A_Inv, 7, 7, -0.7702);

        case 2:
            gsl_matrix_set(A, 8, 8, -1.9964);
            gsl_matrix_set(A, 8, 9, -2.5921);
            gsl_matrix_set(A_Inv, 9, 8, -0.3858);
            gsl_matrix_set(A_Inv, 9, 9, -0.7702);
        case 3:
            gsl_matrix_set(A, 10, 10, -0.9964);
            gsl_matrix_set(A, 10, 11, -2.5921);
            gsl_matrix_set(A_Inv, 11, 10, -0.3858);
            gsl_matrix_set(A_Inv, 11, 11, -0.7702);
        default:
            break;
    }
}



// TODO clearify what index of C_in must be stored in C_op
/*
 * For better understanding look up in Jan's linearized model line 72 to 74
 *          ,_____________________________________________________________,
 *          | 0  C_in1  0  C_in2  0  C_in3  0  C_in4  0  C_in5  0  C_in6  |
 *      C = | 0  C_in7  0  C_in8  0  C_in9  0  C_in10 0  C_in11 0  C_in12 |
 *          | 0  C_in13 0  C_in14 0  C_in15 0  C_in16 0  C_in17 0  C_in18 |
 *          +-------------------------------------------------------------+
 *
 *          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *          !!!!!!!!!!!!!!!!!!!!!!!!! Attention !!!!!!!!!!!!!!!!!!!!!!!!!!!
 *          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *          !!!!!                                                     !!!!!
 *          !!!!!   C_in is very likely a time-dependent vector with  !!!!!
 *          !!!!!   a size of 18. This should not lead to any weird   !!!!!
 *          !!!!!   confusion. So get the data out of the vector      !!!!!
 *          !!!!!                                                     !!!!!
 *          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */
void calculate_Cop(){

    // Setting the values of C_in in the Matrix C_op
    gsl_matrix_set(C_op, 0, 1, gsl_matrix_get(C_in, 0, 0));
    gsl_matrix_set(C_op, 0, 3, gsl_matrix_get(C_in, 1, 0));
    gsl_matrix_set(C_op, 0, 5, gsl_matrix_get(C_in, 2, 0));
    gsl_matrix_set(C_op, 0, 7, gsl_matrix_get(C_in, 3, 0));
    gsl_matrix_set(C_op, 0, 9, gsl_matrix_get(C_in, 4, 0));
    gsl_matrix_set(C_op, 0,11, gsl_matrix_get(C_in, 5, 0));
    gsl_matrix_set(C_op, 1, 1, gsl_matrix_get(C_in, 6, 0));
    gsl_matrix_set(C_op, 1, 3, gsl_matrix_get(C_in, 7, 0));
    gsl_matrix_set(C_op, 1, 5, gsl_matrix_get(C_in, 8, 0));
    gsl_matrix_set(C_op, 1, 7, gsl_matrix_get(C_in, 9, 0));
    gsl_matrix_set(C_op, 1, 9, gsl_matrix_get(C_in,10, 0));
    gsl_matrix_set(C_op, 1,11, gsl_matrix_get(C_in,11, 0));
    gsl_matrix_set(C_op, 2, 1, gsl_matrix_get(C_in,12, 0));
    gsl_matrix_set(C_op, 2, 3, gsl_matrix_get(C_in,13, 0));
    gsl_matrix_set(C_op, 2, 5, gsl_matrix_get(C_in,14, 0));
    gsl_matrix_set(C_op, 2, 7, gsl_matrix_get(C_in,15, 0));
    gsl_matrix_set(C_op, 2, 9, gsl_matrix_get(C_in,16, 0));
    gsl_matrix_set(C_op, 2,11, gsl_matrix_get(C_in,17, 0));


}
/*
 * For better understanding look up in Jan's linearized model line 76 to 78
 *          ,_____________________________________________________________,
 *          | 0  0  0  0  0  0  D_in1  D_in2  D_in3  D_in4  D_in5  D_in6  |
 *      D = | 0  0  0  0  0  0  D_in7  D_in8  D_in9  D_in10 D_in11 D_in12 |
 *          | 0  0  0  0  0  0  D_in13 D_in14 D_in15 D_in16 D_in17 D_in18 |
 *          +-------------------------------------------------------------+
 *
 *          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *          !!!!!!!!!!!!!!!!!!!!!!!!! Attention !!!!!!!!!!!!!!!!!!!!!!!!!!!
 *          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *          !!!!!                                                     !!!!!
 *          !!!!!   D_in is very likely a time-dependent vector with  !!!!!
 *          !!!!!   a size of 18. This should not lead to any weird   !!!!!
 *          !!!!!   confusion. So get the data out of the vector      !!!!!
 *          !!!!!                                                     !!!!!
 *          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */
void calculate_Dop(){

    // Setting the values of D_in in the Matrix D_op
    gsl_matrix_set(D_op, 0, 6, gsl_matrix_get(D_in, 0, 0));
    gsl_matrix_set(D_op, 0, 7, gsl_matrix_get(D_in, 1, 0));
    gsl_matrix_set(D_op, 0, 8, gsl_matrix_get(D_in, 2, 0));
    gsl_matrix_set(D_op, 0, 9, gsl_matrix_get(D_in, 3, 0));
    gsl_matrix_set(D_op, 0,10, gsl_matrix_get(D_in, 4, 0));
    gsl_matrix_set(D_op, 0,11, gsl_matrix_get(D_in, 5, 0));
    gsl_matrix_set(D_op, 1, 6, gsl_matrix_get(D_in, 6, 0));
    gsl_matrix_set(D_op, 1, 7, gsl_matrix_get(D_in, 7, 0));
    gsl_matrix_set(D_op, 1, 8, gsl_matrix_get(D_in, 8, 0));
    gsl_matrix_set(D_op, 1, 9, gsl_matrix_get(D_in, 9, 0));
    gsl_matrix_set(D_op, 1,10, gsl_matrix_get(D_in,10, 0));
    gsl_matrix_set(D_op, 1,11, gsl_matrix_get(D_in,11, 0));
    gsl_matrix_set(D_op, 2, 6, gsl_matrix_get(D_in,12, 0));
    gsl_matrix_set(D_op, 2, 7, gsl_matrix_get(D_in,13, 0));
    gsl_matrix_set(D_op, 2, 8, gsl_matrix_get(D_in,14, 0));
    gsl_matrix_set(D_op, 2, 9, gsl_matrix_get(D_in,15, 0));
    gsl_matrix_set(D_op, 2,10, gsl_matrix_get(D_in,16, 0));
    gsl_matrix_set(D_op, 2,11, gsl_matrix_get(D_in,17, 0));


}

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

