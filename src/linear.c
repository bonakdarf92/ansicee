//
// Created by Farid Bonakdar on 30.06.16.
//

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "linear.h"


gsl_matrix* A;              // Declaration of Matrix A
gsl_matrix* B;              // Declaration of Matrix B
gsl_matrix* A_Inv;          // Declaration of Matrix A inverse
gsl_matrix* KI;             // Declaration of Matrix KI
gsl_matrix* KS;             // Declaration of Matrix KS
gsl_matrix* KP;
gsl_matrix* Ai;             // Declaration of Matrix A iterator
gsl_matrix* Ag;
gsl_matrix* EW_I;           // Declaration of Matrix Eigenvalue for Integrator
gsl_matrix* EW_G;           // Declaration of Matrix Eigenvalue for G??
gsl_matrix* EW_I1;          // Declaration of Matrix Eigenvalue for inner calculations
gsl_matrix* EW_G1;
gsl_matrix* EWI_output;      // Declaration of output Matrix EWI
gsl_matrix* EWG_output;

// TODO check what kind of data C_in and D_in are: Vector or Matrix
gsl_vector* C_in;           // Declaration of input Matrix C
gsl_vector* D_in;           // Declaration of input Matrix D

gsl_matrix* C_op;           // Declaration of operation Matrix C
gsl_matrix* D_op;           // Declaration of operation Matrix D
gsl_matrix* temp1;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temp2;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temp3;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temp4;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temporary;
gsl_matrix* eye;            // Declaration of identity matrix for calculations KP
double a_min = 0.001;       // Declaration and initialization of tuning factor a_min
double A0 = 10;
double b_min = 0.00001;     // Declaration and initialization of tuning factor b_min
double B0 = 10;             // TODO check which value B0 has

gsl_eigen_nonsymm_workspace* workspace;     // Declaration of workspace for eigenvalue calculation
gsl_vector_complex* eigenvalue;             // Declaration of vector for eigenvalues
gsl_vector_complex* eigenvalue2;            // Declaration of vector for eigenvalues
gsl_vector_complex* eigenvalue3;            // Declaration of vector for eigenvalues
gsl_vector_complex* eigenvalue4;            // Declaration of vector for eigenvalues


gsl_vector_view* a;
double delta_a [] = {1.8, 1.3, 1.05};
double delta_b [] = {5, 1.8, 1.2, 1.02};


/*
 * Initialing all matrices and setting up some matrix with zero for easy setup in further code
 * These steps are necessary because the calculations are made every iteration
 */
void initMatrix(){
    A = gsl_matrix_alloc(12,12);
    B = gsl_matrix_alloc(12,12);
    A_Inv = gsl_matrix_alloc(12,12);
    Ai = gsl_matrix_alloc(15,15);
    Ag = gsl_matrix_alloc(15,15);
    KS = gsl_matrix_alloc(12,3);
    KI = gsl_matrix_alloc(12,3);
    KP = gsl_matrix_alloc(12,3);
    C_op = gsl_matrix_alloc(12,3);
    D_op = gsl_matrix_alloc(12,3);
    C_in = gsl_vector_alloc(18);
    D_in = gsl_vector_alloc(18);
    temp1 = gsl_matrix_alloc(12,12);
    temp2 = gsl_matrix_alloc(12,3);
    temp3 = gsl_matrix_alloc(3,3);
    temp4 = gsl_matrix_alloc(3,12);
    temporary = gsl_matrix_alloc(12,3);
    eye = gsl_matrix_alloc(12,12);
    gsl_matrix_set_identity(eye);
    workspace = gsl_eigen_nonsymm_alloc(15);
    eigenvalue = gsl_vector_complex_alloc(15);
    eigenvalue3 = gsl_vector_complex_alloc(15);

    // TODO initialize all matrices


    gsl_matrix_set_zero(A);
    gsl_matrix_set_zero(B);
    gsl_matrix_set_zero(A_Inv);
    gsl_matrix_set_zero(C_op);
    gsl_matrix_set_zero(D_op);

}


/*
 *
 */
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
 * This function is written as a getter
 * It returns the called Matrix
 * 1  --> A
 * 2  --> A_inv
 * 3  --> B
 * 4  --> KS
 * 5  --> KI
 * 6  --> C_in
 * 7  --> C_op
 * 8  --> D_in
 * 9  --> D_op
 * 10 --> Ai
 * 11 --> EW_G
 * 12 --> EW_I
 * 13 --> EWI_end
 * 14 --> KP
 */
gsl_matrix * get_Matrix(size_t n){
    switch (n){
        case 1:
            return A;
        case 2:
            return A_Inv;
        case 3:
            return B;
        case 4:
            return KS;
        case 5:
            return KI;
        //case 6:
          //  return C_in;
        case 7:
            return C_op;
        //case 8:
          //  return D_in;
        case 9:
            return D_op;
        case 10:
            return Ai;
        case 11:
            return EW_G;
        case 12:
            return EW_I;
        case 13:
            EWI_output = EW_I1;
            return EWI_output;
        case 14:
            return KP;
        default:
            return 0;
    }

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
    gsl_matrix_mul_elements(temp2, temp1);
    gsl_matrix_sub(D_op, temp2);
}
// TODO might be a BOTTLENECK !!!


/*
 * This method calculates the matrix -B * Ki (T) and copy all its elements
 * in to the Matrix Ai combined with System matrix A, C and zeros
 *
 *
 *              +---------------------------------------------------------------------------------------+
 *              | A11   A12   A13   A14   A15   A16   A17   A18   A19  A110  A111  A112 | T11  T12  T13 |
 *              | A21   A22   A23   A24   A25   A26   A27   A28   A29  A210  A211  A212 | T21  T22  T23 |
 *              | A31         A33   A34   A35   A36   A37   A38   A39  A310  A311  A312 | T31  T32  T33 |
 *              | A41               A44                                                 | T41  T42  T43 |
 *              | A51                     A55             . . .                 .       | T51  T52  T53 |
 *              | A61         .                 A66                             .       | T61  T62  T63 |
 *              | A71         .                       A77                       .       | T71  T72  T73 |
 *              | A81         .                             A88                         | T81  T82  T83 |
 *              | A91                           .                 A99                   | T91  T92  T93 |
 *              | A101                          .                      A1010            | T101 T102 T103|
 *              | A111                          .                            A1111      | T111 T112 T113|
 *              | A121        . . .                  . . .                        A1212 | T121 T122 T123|
 *              |---------------------------------------------------------------------------------------+
 *              | C11   C12   C13   C14   C15   C16   C17   C18   C19  C110  C111  C112 |  0    0     0 |
 *              | C21   C22   C23   C24   C25   C26   C27   C28   C29  C210  C211  C212 |  0    0     0 |
 *              | C31   C32   C33   C34   C35   C36   C37   C38   C39  C310  C311  C312 |  0    0     0 |
 *              +---------------------------------------------------------------------------------------+
 *
 */
void calculate_Ai(){
    // Scheint zu funktionieren

    gsl_matrix* educt = (gsl_matrix*)getMatrix(5)->data;
    gsl_matrix* product = (gsl_matrix*) getMatrix(3)->data;
    int testVec = gsl_matrix_mul_elements(product, educt);
    printf("%d", testVec);
    gsl_matrix_scale(product, -1);
    temporary = product;

    // Copying the System matrix A (12 x 12) into top left corner of Ai (15 x 15)
    for (size_t i = 0; i < 12; i++) {
        for (size_t j = 0; j < 12; j++) {
            gsl_matrix_set(Ai,i, j, gsl_matrix_get(A, i, j));
        }
    }

    // Copying the temporary stored matrix (-B * Ki) into the top right corner of Ai
    for (size_t k = 0; k < 12 ; k++) {
        for (size_t m = 0; m < 3 ; m++) {
            gsl_matrix_set(Ai, k, 12 + m, gsl_matrix_get(temporary, k, m));
        }
    }

    // Copying the elements of Matrix C_op in bottom left corner of Ai
    for (size_t l = 0; l < 3 ; l++) {
        for (size_t n = 0; n < 12; n++) {
            gsl_matrix_set(Ai, 12 + l, n, gsl_matrix_get(C_op, l, n));
        }
    }

    // Setting zeros in the bottom right corner of Ai
    for (size_t v = 0; v < 3; v++) {
        for (size_t w = 0; w < 3; w++) {
            gsl_matrix_set(Ai, 12+v, 12+w, 0);
        }
    }
}
//TODO Test this method properly


/*
 * This method calculates the eigenvalues of Matrix Ai
 * and stores the data in complex vector eigenvalue
 * The computation is done with help of a computational
 * workspace
 */
void calculate_EWI(){
    EW_I = Ai;
    gsl_eigen_nonsymm(EW_I, eigenvalue, workspace);
}
// TODO Test the method and compare with matlab data


/*
 * This method gets an matrix KS and an integer a0
 * and calculates the matrix KI with the formula taken
 * from Jan's linearized modell line 94
 *
 *      --->  Ki = a0 * (KS' / (KS * KS') )
 *  and stored in the matrix KI
 */
void calculate_KI(gsl_matrix* KS, double a0){

    // temporary storage of matrix KS for later usage
    temp2 = KS;

    // Transpose the matrix KS and save into temp4
    gsl_matrix_transpose_memcpy(temp4,KS);

    // Calculating the product of KS and KS_transpose
    temp3 = KS;
    gsl_matrix_mul_elements(temp3, temp4);

    // Dividing the matrix KS with the Product of KS and KS_transpose
    KI = KS;
    gsl_matrix_div_elements(KI, temp3);

    // Scaling KI with a0
    gsl_matrix_scale(KI, a0);
}
// TODO might be a BOTTLENECK !!!


/*
 * This method changes the matrix indexes for direction change in engine speed
 * For example the direction of movement converts from forward to reverse
 * Then the changes have to be made in System matrix A and A_inverse
 */
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
// TODO check what information is given in n_up down and how the differential observation takes place


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
    gsl_matrix_set(C_op, 0, 1, gsl_vector_get(C_in, 0));
    gsl_matrix_set(C_op, 0, 3, gsl_vector_get(C_in, 1));
    gsl_matrix_set(C_op, 0, 5, gsl_vector_get(C_in, 2));
    gsl_matrix_set(C_op, 0, 7, gsl_vector_get(C_in, 3));
    gsl_matrix_set(C_op, 0, 9, gsl_vector_get(C_in, 4));
    gsl_matrix_set(C_op, 0,11, gsl_vector_get(C_in, 5));
    gsl_matrix_set(C_op, 1, 1, gsl_vector_get(C_in, 6));
    gsl_matrix_set(C_op, 1, 3, gsl_vector_get(C_in, 7));
    gsl_matrix_set(C_op, 1, 5, gsl_vector_get(C_in, 8));
    gsl_matrix_set(C_op, 1, 7, gsl_vector_get(C_in, 9));
    gsl_matrix_set(C_op, 1, 9, gsl_vector_get(C_in,10));
    gsl_matrix_set(C_op, 1,11, gsl_vector_get(C_in,11));
    gsl_matrix_set(C_op, 2, 1, gsl_vector_get(C_in,12));
    gsl_matrix_set(C_op, 2, 3, gsl_vector_get(C_in,13));
    gsl_matrix_set(C_op, 2, 5, gsl_vector_get(C_in,14));
    gsl_matrix_set(C_op, 2, 7, gsl_vector_get(C_in,15));
    gsl_matrix_set(C_op, 2, 9, gsl_vector_get(C_in,16));
    gsl_matrix_set(C_op, 2,11, gsl_vector_get(C_in,17));


}
// TODO clarify what index of C_in must be stored in C_op


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
    gsl_matrix_set(D_op, 0, 6, gsl_vector_get(D_in, 0));
    gsl_matrix_set(D_op, 0, 7, gsl_vector_get(D_in, 1));
    gsl_matrix_set(D_op, 0, 8, gsl_vector_get(D_in, 2));
    gsl_matrix_set(D_op, 0, 9, gsl_vector_get(D_in, 3));
    gsl_matrix_set(D_op, 0,10, gsl_vector_get(D_in, 4));
    gsl_matrix_set(D_op, 0,11, gsl_vector_get(D_in, 5));
    gsl_matrix_set(D_op, 1, 6, gsl_vector_get(D_in, 6));
    gsl_matrix_set(D_op, 1, 7, gsl_vector_get(D_in, 7));
    gsl_matrix_set(D_op, 1, 8, gsl_vector_get(D_in, 8));
    gsl_matrix_set(D_op, 1, 9, gsl_vector_get(D_in, 9));
    gsl_matrix_set(D_op, 1,10, gsl_vector_get(D_in,10));
    gsl_matrix_set(D_op, 1,11, gsl_vector_get(D_in,11));
    gsl_matrix_set(D_op, 2, 6, gsl_vector_get(D_in,12));
    gsl_matrix_set(D_op, 2, 7, gsl_vector_get(D_in,13));
    gsl_matrix_set(D_op, 2, 8, gsl_vector_get(D_in,14));
    gsl_matrix_set(D_op, 2, 9, gsl_vector_get(D_in,15));
    gsl_matrix_set(D_op, 2,10, gsl_vector_get(D_in,16));
    gsl_matrix_set(D_op, 2,11, gsl_vector_get(D_in,17));


}
// TODO clarify what index of D_in must be stored in D_op


/*
 * This method calculates the system matrix Ewi
 * It does the calculations with obtaining the
 * maximum real part of the eigenvalues of the system
 * and adapt it via iteration.
 */
void matrix_Calculator_EWI(){

    // Counter for termination of while loop
    size_t i = 0;

    // Pointer on vector_view containing all real parts of complex vector eigenvalue
    *a = gsl_vector_complex_real(eigenvalue);

    // Storing the structural element vector of vector_view a into gsl_vector real
    gsl_vector* real = &a->vector;

    // Finding the biggest value of vector real
    double max = gsl_vector_max(real);

    /*
     * This while loop calculates all Matrices via iteration of tuning factor a0
     * and the maximum real component.
     * This leads to adaptation of the system matrix in respect to its poles
     *
     */
    while ((max > 0) && (A0 >= a_min) && (i <= I_MAX_A)){
        // If a0 is to big get the half for next iteration
        A0 *= 0.5;

        // Function calls for the iteration
        calculate_KI(KS, A0);
        calculate_Ai();
        calculate_EWI();

        // Increment counter for while loop
        i++;
    }
}
//TODO Test if pointer or address is needed


/*
 * This method do shit
 */
void tune_matrix_EWI(){
    // Initialize matrix EW_I1 with ground matrix EW_I
    EW_I1 = EW_I;
    gsl_eigen_nonsymm(EW_I1,eigenvalue2,workspace);

    // counter for while loop
    size_t j = 0;
    size_t i = 1;

    // The biggest real number of complex eigenvalue vector of EW_I
    const gsl_vector saving = gsl_vector_complex_real(eigenvalue).vector;
    double real1 = gsl_vector_max(&saving);

    // The biggest real number of complex eigenvalue vector of EW_I1 which represents the current state
    const gsl_vector saving2 = gsl_vector_complex_real(eigenvalue2).vector;
    double real2 = gsl_vector_max(&saving2);

    // This for loop iterates over the constant factors of delta_a and compute the while loop
    for (size_t c = 0; c < 3; c++) {

        /*
         * If real1 is negative and real1 of current state is smaller the real2 of previous
         * state and counter is less then 10 then
         */
        while ((real1 <= 0) && (real1 <= real2) && (i+j <= I_MAX_A)){
            A0 *= delta_a[c];       // Adapt Factor A0
            EW_I1 = EW_I;           // Save current Eigenvalues to new one
            calculate_KI(KS,A0);    // Calculates the new matrix KS
            calculate_Ai();         // Calculates the new System matrix Ai
            calculate_EWI();        // Calculates the new eigenvalues
            j++;                    // Increment the counter for while loop
        }

        // Undo the previous multiplication in the while loop for next computation
        A0 /= delta_a[c];

        // Save the new calculated Matrix into the old one
        EW_I = EW_I1;
    }
}
// TODO Very Important!!! Test this method


/*
 * Same method as calculate_KI, instead of using b0 as tuningfactor
 * and the second part calculating the new KP matrix
 */
void calculate_KP(gsl_matrix* KS, double b0){

    // temporary storage of matrix KS for later usage
    temp2 = KS;

    // Transpose the matrix KS and save into temp4
    gsl_matrix_transpose_memcpy(temp4,KS);

    // Calculating the product of KS and KS_transpose
    temp3 = KS;
    gsl_matrix_mul_elements(temp3, temp4);

    // Dividing the matrix KS with the Product of KS and KS_transpose
    KP = KS;
    gsl_matrix_div_elements(KP, temp3);

    // Scaling KI with b0
    gsl_matrix_scale(KP, b0);

    // Saving the current matrix of KP into KP_current for calculations
    gsl_matrix* KP_current = get_Matrix(14);

    // Saving the current subtrahend into sub for calculations
    gsl_matrix_mul_elements(KP_current, D_op);
    gsl_matrix* sub = eye;
    gsl_matrix_sub(sub, KP_current);

    // Calculating the Division of sub and KP_current
    KP = sub;
    gsl_matrix_div_elements(KP, KP_current);
}
// TODO check if its working

/*
 * The structure of this method is nearly the same as calculate_Ai.
 * It computes the System matrix Ag with respect to (A - B * Kp * C) -> [L], (-B * Ki) -> [T], C and zeros
 *
 *              +---------------------------------------------------------------------------------------+
 *              | L11   L12   L13   L14   L15   L16   L17   L18   L19  L110  L111  L112 | T11  T12  T13 |
 *              | L21   L22   L23   L24   L25   L26   L27   L28   L29  L210  L211  L212 | T21  T22  T23 |
 *              | L31         L33   L34   L35   L36   L37   L38   L39  L310  L311  L312 | T31  T32  T33 |
 *              | L41               L44                                                 | T41  T42  T43 |
 *              | L51                     L55             . . .                 .       | T51  T52  T53 |
 *              | L61         .                 L66                             .       | T61  T62  T63 |
 *              | L71         .                       L77                       .       | T71  T72  T73 |
 *              | L81         .                             L88                         | T81  T82  T83 |
 *              | L91                           .                 L99                   | T91  T92  T93 |
 *              | L101                          .                      L1010            | T101 T102 T103|
 *              | L111                          .                            L1111      | T111 T112 T113|
 *              | L121        . . .                  . . .                        L1212 | T121 T122 T123|
 *              |---------------------------------------------------------------------------------------+
 *              | C11   C12   C13   C14   C15   C16   C17   C18   C19  C110  C111  C112 |  0    0     0 |
 *              | C21   C22   C23   C24   C25   C26   C27   C28   C29  C210  C211  C212 |  0    0     0 |
 *              | C31   C32   C33   C34   C35   C36   C37   C38   C39  C310  C311  C312 |  0    0     0 |
 *              +---------------------------------------------------------------------------------------+
 *
 */
void calculate_AG(){

    // Storing the Matrix B into the new Matrix B_current for further calculations
    gsl_matrix* B_current = get_Matrix(3);

    // Copying the matrix B_current in B_current2 for second calculation
    gsl_matrix* B_current2 = B_current;

    // Storing the matrix KI in new matrix KI_current
    gsl_matrix* KI_current = get_Matrix(5);

    // Computing temp_calc1: -->  temp_calc1 = B_current * KI_current
    gsl_matrix* temp_calc1 = B_current;
    gsl_matrix_mul_elements(temp_calc1, KI_current);

    gsl_matrix_scale(temp_calc1, -1);           // Multiply temp_calc1 with -1
    gsl_matrix* A_current = get_Matrix(1);       // Storing the Matrix A in new Matrix A_current
    gsl_matrix* KP_current = get_Matrix(14);     // Storing the Matrix KP in new Matrix KP_current
    gsl_matrix* C_current = get_Matrix(7);       // Storing the Matrix C_current in new Matrix C_current

    /*
     * Computing the next three steps can be summarized into this picture:
     *
     *      --> temp_calc2 = KP_current * C_current
     *      --> temp_calc3 = B_current2 * temp_calc2 = B_current * KP_current * C_current
     *      --> temp_calc4 = A_current - temp_calc3 = ...
     *
     *  These steps are necessary to obtain the big Matrix [L]
     *      --> [L] = A - B * KP * C
     *
     */
    gsl_matrix* temp_calc2 = KP_current;
    gsl_matrix_mul_elements(temp_calc2, C_current);
    gsl_matrix* temp_calc3 = B_current2;
    gsl_matrix_mul_elements(B_current2, temp_calc2);
    gsl_matrix* temp_calc4 = A_current;
    gsl_matrix_sub(A_current, temp_calc3);
    // TODO check it --> will probably not work

    // Copying the System matrix L (12 x 12) into top left corner of AG (15 x 15)
    for (size_t i = 0; i < 12; i++) {
        for (size_t j = 0; j < 12; j++) {
            gsl_matrix_set(Ag,i, j, gsl_matrix_get(temp_calc4, i, j));
        }
    }

    // Copying the temporary stored matrix (-B * Ki) into the top right corner of Ag
    for (size_t k = 0; k < 12 ; k++) {
        for (size_t m = 0; m < 3 ; m++) {
            gsl_matrix_set(Ag, k, 12 + m, gsl_matrix_get(temp_calc1, k, m));
        }
    }

    // Copying the elements of Matrix C_op in bottom left corner of Ag
    for (size_t l = 0; l < 3 ; l++) {
        for (size_t n = 0; n < 12; n++) {
            gsl_matrix_set(Ag, 12 + l, n, gsl_matrix_get(C_op, l, n));
        }
    }

    // Setting zeros in the bottom right corner of Ag
    for (size_t v = 0; v < 3; v++) {
        for (size_t w = 0; w < 3; w++) {
            gsl_matrix_set(Ag, 12+v, 12+w, 0);
        }
    }
}
// TODO Standard check


/*
 * Comment
 */
void calculate_EWG(){
    EW_G = Ag;
    gsl_eigen_nonsymm(EW_G, eigenvalue3, workspace);
}
// TODO check this methd, pretty same like calculate_EWI


/*
 * Comment
 */
void matrix_Calculator_EWG(){
    // Counter for termination of while loop
    size_t i = 0;

    // Pointer on vector_view containing all real parts of complex vector eigenvalue3
    *a = gsl_vector_complex_real(eigenvalue3);

    // Storing the structural element vector of vector_view a into gsl_vector real
    gsl_vector* real = &a->vector;

    // Finding the biggest value of vector real
    double max = gsl_vector_max(real);

    /*
     * This while loop calculates all Matrices via iteration of tuning factor a0
     * and the maximum real component.
     * This leads to adaptation of the system matrix in respect to its poles
     *
     */
    while ((max > 0) && (B0 >= b_min) && (i <= I_MAX_B)){
        // If a0 is to big get the half for next iteration
        B0 *= 0.5;

        // Function calls for the iteration
        calculate_KP(KS, B0);
        calculate_AG();
        calculate_EWG();

        // Increment counter for while loop
        i++;
    }
}
// TODO check it


/*
 * Comment
 */
void tune_matrix_EWG(){

    // Initialize matrix EW_I1 with ground matrix EW_I
    EW_G1 = EW_G;
    gsl_eigen_nonsymm(EW_G1,eigenvalue4,workspace);

    // counter for while loop
    size_t j = 0;
    size_t i = 1;

    // The biggest real number of complex eigenvalue vector of EW_G
    const gsl_vector saving = gsl_vector_complex_real(eigenvalue3).vector;
    double real1 = gsl_vector_max(&saving);

    // The biggest real number of complex eigenvalue vector of EW_G1 which represents the current state
    const gsl_vector saving2 = gsl_vector_complex_real(eigenvalue4).vector;
    double real2 = gsl_vector_max(&saving2);

    // This for loop iterates over the constant factors of delta_a and compute the while loop
    for (size_t c = 0; c < 4; c++) {

        /*
         * If real1 is negative and real1 of current state is smaller the real2 of previous
         * state and counter is less then 10 then
         */
        while ((real1 <= 0) && (real1 <= real2) && (i+j <= I_MAX_B)){
            B0 *= delta_b[c];       // Adapt Factor A0
            EW_G1 = EW_G;           // Save current Eigenvalues to new one
            calculate_KP(KS,B0);    // Calculates the new matrix KS
            calculate_AG();         // Calculates the new System matrix Ai
            calculate_EWG();        // Calculates the new eigenvalues
            j++;                    // Increment the counter for while loop
        }

        // Undo the previous multiplication in the while loop for next computation
        B0 /= delta_b[c];

        // Save the new calculated Matrix into the old one
        EW_G = EW_G1;
    }
}
// TODO check it


/*
 * Comment
 */
void tune_KP(){
    double b_end = B0;
    EWG_output = EW_G1;
    calculate_KP(KS, b_end);

}
// TODO check