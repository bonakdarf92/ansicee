//
// Created by Farid Bonakdar on 30.06.16.
//

#include <InitTest.h>
#include <calcMatrix_A_B.h>
#include <gsl/gsl_matrix.h>
#include "linear.h"
#include <gsl/gsl_errno.h>


gsl_matrix* ASystem;        // Declaration of Matrix A
gsl_matrix* BSystem;        // Declaration of Matrix B
gsl_matrix* A_Inv;          // Declaration of Matrix A inverse
gsl_matrix* KI;             // Declaration of Matrix KI
gsl_matrix* KS;             // Declaration of Matrix KS
gsl_matrix* KP;
gsl_matrix* Ai;             // Declaration of Matrix A iterator
gsl_matrix* Ag;
gsl_matrix_complex* EW_I;           // Declaration of Matrix Eigenvalue for Integrator
gsl_matrix* EW_G;           // Declaration of Matrix Eigenvalue for G??
gsl_matrix* EW_I1;          // Declaration of Matrix Eigenvalue for inner calculations
gsl_matrix* EW_G1;
gsl_matrix* EWI_output;     // Declaration of output Matrix EWI
gsl_matrix* EWG_output;

gsl_vector* C_in;           // Declaration of input Matrix C
gsl_vector* D_in;           // Declaration of input Matrix D

gsl_matrix* C_op;           // Declaration of operation Matrix C
gsl_matrix* D_op;           // Declaration of operation Matrix D
gsl_matrix* temp1;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temp2;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temp3;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temp4;          // Declaration of temporary Matrix for Calculations
gsl_matrix* temp5;
gsl_matrix* temp6;
gsl_matrix* temporary;
gsl_matrix* Inverse;
gsl_matrix* eye;            // Declaration of identity matrix for calculations KP
double a_min = 0.001;       // Declaration and initialization of tuning factor a_min
double A0 = 10.0;/**/
double b_min = 0.00001;     // Declaration and initialization of tuning factor b_min
double B0 = 10.0;             // TODO check which value B0 has

gsl_eigen_nonsymm_workspace* workspace;     // Declaration of workspace for eigenvalue calculation
gsl_vector_complex* eigenvalue;             // Declaration of vector for eigenvalues
gsl_vector_complex* eigenvalue2;            // Declaration of vector for eigenvalues
gsl_vector_complex* eigenvalue3;            // Declaration of vector for eigenvalues
gsl_vector_complex* eigenvalue4;            // Declaration of vector for eigenvalues
gsl_matrix* dmd;
gsl_matrix* threetimes3;
gsl_matrix* linalgtest;
gsl_matrix* KI2;
gsl_matrix* KI3;
gsl_matrix* KI4;

gsl_vector_view a;
double delta_a [] = {1.8, 1.3, 1.05};
double delta_b [] = {5, 1.8, 1.2, 1.02};


/*
 * Initialing all matrices and setting up some matrix with zero for easy setup in further code
 * These steps are necessary because the calculations are made every iteration
 */
void initMatrix(void){
    ASystem = gsl_matrix_alloc(12,12);
    BSystem = gsl_matrix_alloc(12,12);
    A_Inv = gsl_matrix_alloc(12,12);
    Ai = gsl_matrix_alloc(15,15);
    Ag = gsl_matrix_alloc(15,15);
    KS = gsl_matrix_alloc(3,12);
    KI = gsl_matrix_alloc(12,3);
    KP = gsl_matrix_alloc(12,3);
    C_op = gsl_matrix_alloc(3,12);
    D_op = gsl_matrix_alloc(3,12);
    C_in = gsl_vector_alloc(18);
    D_in = gsl_vector_alloc(18);
    temp1 = gsl_matrix_alloc(12,12);
    //temp2 = gsl_matrix_alloc(12,3);
    temp2 = gsl_matrix_alloc(3,12);
    temp3 = gsl_matrix_alloc(3,3);
    temp4 = gsl_matrix_alloc(3,12);
    temp5 = gsl_matrix_alloc(12,3);
    temp6 = gsl_matrix_alloc(12, 12);
    temporary = gsl_matrix_alloc(12,3);
    eye = gsl_matrix_alloc(12,12);
    gsl_matrix_set_identity(eye);
    workspace = gsl_eigen_nonsymm_alloc(15);
    eigenvalue = gsl_vector_complex_alloc(15);
    eigenvalue2 = gsl_vector_complex_alloc(15);
    eigenvalue3 = gsl_vector_complex_alloc(15);
    eigenvalue4 = gsl_vector_complex_alloc(15);
    dmd = gsl_matrix_alloc(3,3);
    Inverse = gsl_matrix_alloc(3, 3);
    EW_I = gsl_matrix_complex_alloc(15, 15);
    EW_I1 = gsl_matrix_alloc(15, 15);
    KI2 = gsl_matrix_alloc(12, 3);
    KI3 = gsl_matrix_alloc(12, 3);
    KI4 = gsl_matrix_alloc(12, 3);
    threetimes3 = gsl_matrix_alloc(3, 3);
    linalgtest = gsl_matrix_alloc(3, 3);



    // TODO initialize all matrices
    gsl_matrix_set_zero(ASystem);
    gsl_matrix_set_zero(BSystem);
    gsl_matrix_set_zero(A_Inv);
    gsl_matrix_set_zero(C_op);
    gsl_matrix_set_zero(D_op);

}

/*
 * This Method do the presetting of the matrices
 * The numbers and the shape is taken from Jan's linearized Model
 * For further details look in the comments or the source code
 * in Matlab Simulink
 */
void matrixPresetting(void){
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
    gsl_matrix_set(ASystem, 0, 0, -7.2384);           // Setting value in Matrix A_1_1
    gsl_matrix_set(ASystem, 0, 1, -38.9376);          // Setting value in Matrix A_1_2
    gsl_matrix_set(ASystem, 1, 0, 1);                 // Setting value in Matrix A_2_1
    gsl_matrix_set(ASystem, 2, 2, -7.2384);           // Setting value in Matrix A_3_3
    gsl_matrix_set(ASystem, 2, 3, -38.9376);          // Setting value in Matrix A_3_4
    gsl_matrix_set(ASystem, 3, 2, 1);                 // Setting value in Matrix A_4_3
    gsl_matrix_set(ASystem, 4, 4, -7.2384);           // Setting value in Matrix A_5_5
    gsl_matrix_set(ASystem, 4, 5, -38.9376);          // Setting value in Matrix A_5_6
    gsl_matrix_set(ASystem, 5, 4, 1);                 // Setting value in Matrix A_6_5
    gsl_matrix_set(ASystem, 6, 6, -7.0356);           // Setting value in Matrix A_7_7
    gsl_matrix_set(ASystem, 6, 7, -183.0609);         // Setting value in Matrix A_7_8
    gsl_matrix_set(ASystem, 7, 6, 1);                 // Setting value in Matrix A_8_7
    gsl_matrix_set(ASystem, 8, 8, -7.0356);           // Setting value in Matrix A_9_9
    gsl_matrix_set(ASystem, 8, 9, -183.0609);         // Setting value in Matrix A_9_10
    gsl_matrix_set(ASystem, 9, 8, 1);                 // Setting value in Matrix A_10_9
    gsl_matrix_set(ASystem, 10, 10, -7.0356);         // Setting value in Matrix A_11_11
    gsl_matrix_set(ASystem, 10, 11, -183.0609);       // Setting value in Matrix A_11_12
    gsl_matrix_set(ASystem, 11, 10, 1);               // Setting value in Matrix A_12_11

    //TODO T_n1, T_n2 und T_n3
    // Setting all the indexes for Matrix B
    gsl_matrix_set(BSystem, 0, 0, (KONSTANTE / (T_DELTA*T_DELTA)) );
    gsl_matrix_set(BSystem, 2, 1, (KONSTANTE / (T_DELTA*T_DELTA)) );
    gsl_matrix_set(BSystem, 4, 2, (KONSTANTE / (T_DELTA*T_DELTA)));
    //gsl_matrix_set(BSystem, 6, 3, (KONSTANTE / (T_N_CONSTUP*T_N_CONSTDN)));     // TODO T_n_1
    //gsl_matrix_set(BSystem, 8, 4, (KONSTANTE / (T_N_CONSTDN*T_N_CONSTDN));     // TODO T_n_2
    //gsl_matrix_set(BSystem, 10, 5, (KONSTANTE / (T_N_CONSTDN*T_N_CONSTDN)));    // TODO T_n_3


    // Setting all the indexes for Matrix A inverse
    gsl_matrix_set(A_Inv, 0, 1, 1);             // Setting value in Matrix A_inv 1_2
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
        case 1:/**/
            return ASystem;
        case 2:
            return A_Inv;
        case 3:
            return BSystem;
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
        //case 12:
            //return EW_I;
        case 13:
            EWI_output = EW_I1;
            return EWI_output;
        case 14:
            return KP;
        case 15:
            return temp4;
        case 16:
            return Ag;
        case 17:
            return Inverse;
        default:
            return 0;
    }

}

/*
 * This Method calls the function getMatrix from horizontal model and
 * store its returns in Matrix C_in and D_in.
 */
void getInputParameter(void){
    const gsl_vector* save = getMatrix(1);
    const gsl_vector* save2 = getMatrix(2);
    gsl_vector_memcpy(C_in, save);
    gsl_vector_memcpy(D_in, save2);
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
void calculate_KS(void){
    //KS = D_op;
    //gsl_matrix_memcpy(KS, D_op);
    //temp1 = A_Inv;
    //gsl_matrix_memcpy(temp1, A_Inv);
    //temp2 = C_op;
    //gsl_matrix_memcpy(temp2, C_op);
    //gsl_matrix_mul_elements(temp2, A_Inv);
    //gsl_matrix_mul_elements(temp2, BSystem);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, C_op, A_Inv, 0.0, temp2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp2, BSystem, 0.0, temp4);

    //gsl_matrix_sub(temp2, D_op);
    gsl_matrix_memcpy(KS, D_op);
    //gsl_matrix_scale(KS, 1.0);
    gsl_matrix_sub(KS, temp4);
    //printf("Cop ");
    //printer(C_op, NULL);
    //printf("Ainv ");
    //printer(A_Inv, NULL);
    //printf("Berech. prod1 ");
    //printer(temp4, NULL);
    //printf("Korrekter prod1 ");
    //printer(savingMatrix(32), NULL);
    //printf("Berech. prod2 ");
    //printer(temp4, NULL);
    //printf("Korrekter prod2 ");
    //printer(savingMatrix(33), NULL);
    //printf("Berehnet: ");
    //printer(KS, NULL);

}
// TODO Abweichungen an manchen Stellen bis 2 Nachkommastellen !!!

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
void calculate_Ai(void){
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, BSystem, KI, 0.0, temp5);

    // Copying the System matrix A (12 x 12) into top left corner of Ai (15 x 15)
    for (size_t i = 0; i < 12; i++) {
        for (size_t j = 0; j < 12; j++) {
            gsl_matrix_set(Ai,i, j, gsl_matrix_get(ASystem, i, j));
        }
    }

    // Copying the temporary stored matrix (-B * Ki) into the top right corner of Ai
    for (size_t k = 0; k < 12 ; k++) {
        for (size_t m = 0; m < 3 ; m++) {
            if ( gsl_isnan(gsl_matrix_get(temp5, k, m)) || gsl_isinf(gsl_matrix_get(temp5, k, m))) {
                //gsl_matrix_set(Ai, k, 12 + m, gsl_matrix_get(temp5, k, m));
                gsl_matrix_set(Ai, k, 12 + m, 0);
            } else
                //gsl_matrix_set(Ai, k, 12 + m, 0);
                gsl_matrix_set(Ai, k, 12 + m, gsl_matrix_get(temp5, k, m));
        }
    }

    // Copying the elements of Matrix C_op in bottom left corner of Ai
    for (size_t l = 0; l < 3 ; l++) {
        for (size_t n = 0; n < 12; n++) {
            if ( gsl_isnan(gsl_matrix_get(C_op, l, n)) || gsl_isinf(gsl_matrix_get(C_op, l, n)) ) {
                //gsl_matrix_set(Ai, 12 + l, n, gsl_matrix_get(C_op, l, n));
                gsl_matrix_set(Ai, 12 + l, n, 0);
            } else
                //gsl_matrix_set(Ai, 12 + l, n, 0);
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
//TODO Funktioniert jedoch mit Abweichungen aus Ki

/*
 * This method calculates the eigenvalues of Matrix Ai
 * and stores the data in complex vector eigenvalue
 * The computation is done with help of a computational
 * workspace
 */
void calculate_EWI(void){
    // Initialize the complex Matrix view EWI with the values of Ai
    gsl_matrix_complex_view EWI_view = gsl_matrix_complex_view_array(Ai->data, Ai->size1, Ai->size2);

    // Copy the elements of Matrix View into complex Matrix EWI
    gsl_matrix_complex_memcpy(EW_I, &EWI_view.matrix);
    //gsl_eigen_nonsymm(EW_I, eigenvalue, workspace);

    // Define Settings for zgeev (Complex GEneral EigenValues)
    // N = No calculation V = do computation
    char Eigenvectors = 'N';

    // Size of rows/cols
    int n = (int) EW_I->size1;

    // Eigenvalues
    complex double w;

    // Work for auxilary calculation
    complex double work;

    // Size of work
    int lworker = 30;

    // Work
    double rworker;

    // Status if succesful or aborted
    int info;

    // Initalizing complex Matrix for calculations in zgeev
    complex double * A = (complex double *) EW_I->data;

    // Calculate eigenvalues without eigenvectors
    zgeev_(&Eigenvectors, &Eigenvectors, &n, A, &n, &w, NULL, &n, NULL, &n, &work, &lworker, &rworker, &info);
}
// TODO Test the method and compare with matlab data

/*
 * This method gets an matrix KS and an integer a0
 * and calculates the matrix KI with the formula taken
 * from Jan's linearized model line 94
 *
 *      --->  Ki = a0 * (KS' / (KS * KS') )
 *                      |___|   |________|
 *                      temp5     temp 3
 *
 *
 *  and stored in the matrix KI
 */
void calculate_KI(gsl_matrix* KS, double a0){
/*
    // temporary storage of matrix KS for later usage
    //size_t i;
    //gsl_permutation* p1 = gsl_permutation_alloc(3);
    //gsl_permutation* p2 = gsl_permutation_alloc(3);
    //gsl_permutation* p3 = gsl_permutation_alloc(3);

    // Transpose the matrix KS and save into temp5
    gsl_matrix_transpose_memcpy(temp5, KS);

    // Declaration and Initialization Vector tau for QR Decomposition
    gsl_vector* tau = gsl_vector_alloc(3);

    // Declaration and Initialization of Vector b for linear solution
    gsl_vector* b = gsl_vector_alloc(3);

    // Declaration and Initialization of Vector x for linear solution
    gsl_vector* x = gsl_vector_alloc(3);

    // Declaration and Initialization of Vector x for linear solution
    gsl_vector* norm = gsl_vector_alloc(3);

    gsl_permutation* permutation = gsl_permutation_alloc(3);

    gsl_matrix* Q = gsl_matrix_alloc(3, 3);

    gsl_matrix* R = gsl_matrix_alloc(3, 3);

    int signum;

    //int sig1;
    //int sig2;
    //int sig3;
    //gsl_vector* tau1 = gsl_vector_alloc(3);
    //gsl_vector* tau2 = gsl_vector_alloc(3);
    //gsl_vector* norm = gsl_vector_alloc(3);
    //gsl_vector* b = gsl_vector_alloc(3);
    //gsl_vector* x1 = gsl_vector_alloc(3);
    //gsl_vector* x2 = gsl_vector_alloc(3);
    //gsl_vector* x3 = gsl_vector_alloc(3);
    //gsl_matrix* Q = gsl_matrix_alloc(3, 3);
    //gsl_matrix* Btest = gsl_matrix_calloc(12, 3);

    //gsl_matrix* R = gsl_matrix_alloc(3, 3);
*/
    //printer(KS, NULL);
    // Build matrix-matrix product KS * KS' and save it in temp3
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, KS, KS, 0.0, temp3);
    //gsl_matrix_memcpy(temp5, KS);
    //printer(KS, NULL);
    // Transpose matrix temp3
    //gsl_matrix_transpose(temp3);
    //gsl_matrix_memcpy(dmd, temp3);
    //gsl_matrix_memcpy(Inverse, temp3);
    gsl_matrix_memcpy(threetimes3, temp3);
    //gsl_matrix_transpose_memcpy(Btest, KS);

    size_t i;
    size_t j;

    //gsl_linalg_LU_decomp(temp3, p1, &sig1);

    gsl_error_handler_t * handler = gsl_set_error_handler_off();
/*
    //gsl_linalg_LU_decomp(temp3, p1, &sig1);
    //gsl_linalg_LU_decomp(threetimes3, p3, &sig3);
    //gsl_linalg_QR_decomp(dmd, tau1);
    //gsl_linalg_QRPT_decomp(Inverse, tau2, p2, &sig2, norm);
    //gsl_linalg_QRPT_decomp2(Inverse, Q, R, tau2, p2, &sig2, norm);

    //gsl_linalg_LU_solve(temp3, p1, testB, x1);
    //gsl_linalg_QR_solve(dmd, tau1, testB, x2);
    //gsl_linalg_QRPT_solve(Inverse, tau2, p2, testB, x3);


    // Calculating the QR Decomposition of (KS * KS')' and save data in temp3 and tau --> A = QR
    //gsl_linalg_QR_decomp(temp3, tau);
    //gsl_linalg_QRPT_decomp(temp3, tau, permutation, &signum, norm);
    //gsl_linalg_QRPT_decomp2(temp3, Q, R, tau, permutation, &signum, norm);
    //printer(NULL, tau);
    //printf("QR Decomp");
    //printer(NULL, x2);
    //printf("LU decomp ");
    //printer(NULL, x1);
    //printf("QRP' Decomp");
    //printer(NULL, x3);

    //printf("QR Decomp");
    //printer(dmd, NULL);
    //printf("LU decomp ");
    //printer(temp3, NULL);
    //printf("QRP' Decomp");
    //printer(R, NULL);
    */
    /*
    double tauLap;
    double workLap[36];
    int lworkLap = 3;
    //double matrix[3, 3] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    int info;
    int m = 12;
    int n = 3;
    int lda = 3;
    char trans = 'T';
    int nrhs = 12;
    int lwork = 36;

    char uplo = 'U';
    int ipiv[3];
*/

    //double a = threetimes3->data;  // 3 Reihen, 3 Spalten
    //double b = KS->data;           // 3 Reihen, 12 Spalten
    //gsl_matrix_view view1 = gsl_matrix_view_array(threetimes3->data, 3, 3);
    //double a[25] = {6.80, -2.11, 5.66, 5.97, 8.23, -6.05,-3.30, 5.36,-4.44, 1.08, -0.45, 2.58,-2.70, 0.27, 9.04,
    //                 8.32, 2.71, 4.35,-7.17, 2.14, -9.67,-5.14,-7.26, 6.08,-6.87};
    double a3 [9] = {19.0, 2.0, 3.0, 2.0, 4.0, 5.0, 3.0, 5.0, 20.0};
    double b4 [36] = {397.472512951577, -6931.90068901920, -17104.4041448888, 690.259868205320, -12038.0974786996
                ,-29703.8999328737, 387.899070390892, -6764.94033095593, -16692.4309259111, 9740.85235817155 ,-169879.976534301
                 ,-419177.351944731, 875.929366966719, -15276.1642240868, -37693.8012234337, 797.006091333060, -13899.7462557548
            ,-34297.5019602415, -4.13083627588057, 72.0415774524024, 177.761960429290, 1.31554631783812, -22.9430617965015
          ,-56.6118037308531, 0.555002936395560, -9.67922337230086, -23.8833987669135, 1698.80980988175, -29627.1578735975,
                      -73104.7521695148, 1113.07836687159, -19412.0308872516, -47899.0159357860, -841.344182865860, 14673.0003481289,
                      36205.4996503421};
    double a4[9] = {temp3->data[0], temp3->data[1], temp3->data[2], temp3->data[3], temp3->data[4], temp3->data[5],
                    temp3->data[6], temp3->data[7], temp3->data[8]};
    double b5[36];
    double a5[9];
    for (i = 0; i < 36; i++) {
        b5[i] = KS->data[i];
    }
    for (i = 0;  i<9 ; i++) {
        a5[i] = temp3->data[i];
    }
    double b3 [9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    //double *a2 = threetimes3->data;
    //double b[15] = {4.02, 6.19, -8.22, -7.57, -3.03,
    //                -1.56, 4.00,-8.67,  1.75,  2.86,
    //                9.81, -4.09, -4.57, -8.61, 8.99};
    //double *b2 = KS->data;
    int m = 3;
    int n = 3;
    int nrhs = 12;
    int lda = 3;
    int ldb = 3;
    int info;
    int ipiv;
    double xPtr[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int ldx = 3;
    double *workPtr;
    //float swork [18];
    int lwork = 10;
    int iter;
    char trans = 'N';
    double wkopt;
    double* worker;
    int jpvt;
    double tau[3];
    double work[1920];
    //printf("Vorher");
    //printer(threetimes3, NULL);
    //dgeqrf_(&m, &n, threetimes3->data, &lda, &tauLap, &workLap, &lworkLap, &info);
    //printer(linalgtest, NULL);
    //printer(KS, NULL);
    //printf("Nachher");
    //printer(threetimes3, NULL);

    gsl_matrix* testB = gsl_matrix_alloc(3, 12);
    gsl_matrix_memcpy(testB, KS);
    //dgels_(&trans, &m, &n, &nrhs, a5, &lda, b5, &ldb, &wkopt, &lwork, &info);
    //lwork = (int)wkopt;
    //worker = (double*) malloc(lwork * sizeof(double));
    //dgels_(&trans, &m, &n, &nrhs, a5, &lda, b5, &ldb, &wkopt, &lwork, &info);
    //info = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', m, n, nrhs, temp3->data, lda, threetimes3->data, ldb);
    //dsysv_(&uplo, &n, &nrhs, a3, &lda, &ipiv, b3, &ldb, &workPtr, &lwork, &info);
    //dgelsy_(&m, &n, &nrhs, threetimes3->data, &lda, testB->data, &ldb, jpvtPtr, &rcond, &rank, &work, &lwork, &info);
    //dgesv_(&n, &nrhs, a5, &lda, ipiv, b5, &ldb, &info);
    //dsgesv_(&n, &nrhs, a3, &lda, ipiv, b3, &ldb, xPtr, &ldx, workPtr, &swork, &iter, &info);
    //info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, temp3->data, lda, &ipiv, threetimes3->data, ldb);

    printer(temp3, NULL);

    gsl_matrix* output = gsl_matrix_alloc(3, 12);
    //for (j = 0; j < 36; j++) {
    //    output->data[j] = b5[j];
    //}
    output->data = temp3->data;
    gsl_matrix* output2 = gsl_matrix_alloc(3, 3);
    //gsl_matrix_transpose_memcpy(output2, output);
    printer(temp3, NULL);
    printf(" Info %d",info);
/*
    for (i = 0; i < KS->size2; i++) {       // KS.size2 = 12
        gsl_matrix_get_col(b, KS, i);
        gsl_linalg_LU_solve(temp3, p1, b, x1);
        gsl_linalg_QR_solve(dmd, tau1, b, x2);
        gsl_linalg_QRPT_solve(Inverse, tau2, p2, b, x3);

        // Declaration and Initialization Vector tau for QR Decomposition
        //gsl_vector* tau = gsl_vector_alloc(3);

        // Declaration and Initialization of Vector b for linear solution
        //gsl_vector* b = gsl_vector_alloc(3);

        // Declaration and Initialization of Vector x for linear solution
        //gsl_vector* x = gsl_vector_alloc(3);

        // Declaration and Initialization of Vector x for linear solution
        //gsl_vector* norm = gsl_vector_alloc(3);

        //gsl_permutation* permutation = gsl_permutation_alloc(3);

        //gsl_matrix* Q = gsl_matrix_alloc(3, 3);

        //gsl_matrix* R = gsl_matrix_alloc(3, 3);

        //int signum;

        //gsl_linalg_QRPT_decomp2(temp3, Q, R, tau, permutation, &signum, norm);
        //gsl_matrix_get_col(b, KS, i);
        //printf("\nb ");
        //printer(NULL, b);

        //printf("Matrix R ");
        //printer(R, NULL);
        //printf("Matrix Q ");
        //printer(Q, NULL);
        //gsl_linalg_QR_QRsolve(Q, R, b, x);
        gsl_matrix_set_row(KI, i, x1);
        gsl_matrix_set_row(KI2, i, x2);
        gsl_matrix_set_row(KI3, i, x3);
        //printf("x ");
        //printer(NULL, x);
        //printf("QR Decomp");
        //printer(NULL, x2);
        //printf("LU decomp ");
        //printer(NULL, x1);
        //printf("QRP' Decomp");
        //printer(NULL, x3);
     }
*/
    //printf("QR Decomp");
    //printer(KI2, NULL);
    //printf("LU decomp ");
    //printer(KI, NULL);
    //printf("QRP' Decomp");
    //printer(KI3, NULL);
    //printf("LU Inverse ");
    //printer(KI4, NULL);

    // Dividing the matrix KS with the Product of KS and KS_transpose
    // TODO Berechnung überprüfen Werte Stimmen nicht überein

    gsl_set_error_handler(handler);

    // Scaling KI with a0
    //gsl_matrix_scale(KI, a0);

}
// TODO große Abweichungen in den Matrixeinträgen

/*
 * This method changes the matrix indexes for direction change in engine speed
 * For example the direction of movement converts from forward to reverse
 * Then the changes have to be made in System matrix A and A_inverse
 */
void changing_engineSpeed(gsl_matrix* n_updn, size_t zaehler){

    double updn1 = gsl_matrix_get(n_updn, zaehler, 0);
    double updn2 = gsl_matrix_get(n_updn, zaehler, 1);
    double updn3 = gsl_matrix_get(n_updn, zaehler, 2);
    if (updn1 == 0) {
        gsl_matrix_set(ASystem, 6, 6, -1.9964);
        gsl_matrix_set(ASystem, 6, 7, -2.5921);
        gsl_matrix_set(A_Inv, 7, 6, -0.3858);
        gsl_matrix_set(A_Inv, 7, 7, -0.7702);
    }
    if (updn2 == 0) {
        gsl_matrix_set(ASystem, 8, 8, -1.9964);
        gsl_matrix_set(ASystem, 8, 9, -2.5921);
        gsl_matrix_set(A_Inv, 9, 8, -0.3858);
        gsl_matrix_set(A_Inv, 9, 9, -0.7702);
    }
    if (updn3 == 0) {
        gsl_matrix_set(ASystem, 10, 10, -1.9964);
        gsl_matrix_set(ASystem, 10, 11, -2.5921);
        gsl_matrix_set(A_Inv, 11, 10, -0.3858);
        gsl_matrix_set(A_Inv, 11, 11, -0.7702);
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
void calculate_Cop(void){

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

#if DEBUGGER == 1
    printf("\n\n c Matrix Ueb. :");
    printer(NULL, returnCinDin(0));

    printer(C_op, NULL);
#endif
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
void calculate_Dop(void){

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


#if DEBUGGER == 1
    printf("\n\n d Matrix Ueb. :");
    printer(NULL, returnCinDin(1));

    printer(get_Matrix(9), NULL);
#endif
}

/*
 * This method calculates the system matrix Ewi
 * It does the calculations with obtaining the
 * maximum real part of the eigenvalues of the system
 * and adapt it via iteration.
 */
void matrix_Calculator_EWI(void){

    // Counter for termination of while loop
    size_t i = 0;
    double A0temp = A0;

    // Pointer on vector_view containing all real parts of complex vector eigenvalue
    a = gsl_vector_complex_real(eigenvalue);

    // Storing the structural element vector of vector_view a into gsl_vector real
    //gsl_vector* real = a;//->vector;

    // Finding the biggest value of vector real
    double max = gsl_vector_max(&a.vector);
    //printf(" Wert %f\n ", max);

    /*
     * This while loop calculates all Matrices via iteration of tuning factor a0
     * and the maximum real component.
     * This leads to adaptation of the system matrix in respect to its poles
     *
     */
    while ((max > 0) && (A0temp >= a_min) && (i <= I_MAX_A)){
        // If a0 is to big get the half for next iteration
        A0temp *= 0.5;

        // Function calls for the iteration
        calculate_KI(KS, A0temp);
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
void tune_matrix_EWI(void){
    // Initialize matrix EW_I1 with ground matrix EW_I
    //EW_I1 = EW_I;
    gsl_matrix_memcpy(EW_I1, EW_I);
    gsl_eigen_nonsymm(EW_I1,eigenvalue2,workspace);

    // counter for while loop
    size_t j = 0;
    size_t i = 1;
    double A0temp = 10;

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
            A0temp *= delta_a[c];       // Adapt Factor A0
            EW_I1 = EW_I;           // Save current Eigenvalues to new one
            calculate_KI(KS,A0temp);    // Calculates the new matrix KS
            calculate_Ai();         // Calculates the new System matrix Ai
            calculate_EWI();        // Calculates the new eigenvalues
            j++;                    // Increment the counter for while loop
        }

        // Undo the previous multiplication in the while loop for next computation
        A0temp /= delta_a[c];

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
    size_t i, j;

    // Transpose the matrix KS and save into temp5
    gsl_matrix_transpose_memcpy(temp5, KS);

    // Declaration and Initialization Vector tau for QR Decomposition
    gsl_vector* tau = gsl_vector_alloc(3);

    // Declaration and Initialization of Vector b for linear solution
    gsl_vector* b = gsl_vector_alloc(3);

    // Declaration and Initialization of Vector x for linear solution
    gsl_vector* x = gsl_vector_alloc(3);

    // Declaration and Initialization of Vector x for linear solution
    gsl_vector* norm = gsl_vector_alloc(3);

    gsl_permutation* permutation = gsl_permutation_alloc(3);

    int signum;

    // Transpose the matrix KS and save into temp4
    gsl_matrix_transpose_memcpy(temp5,KS);

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, KS, KS, 0.0, temp3);

    gsl_matrix_transpose(temp3);

    // Dividing the matrix KS with the Product of KS and KS_transpose

    gsl_linalg_QRPT_decomp(temp3, tau, permutation, &signum, norm);
    //printer(NULL, tau);

    for (i = 0; i < KS->size2; i++) {
        gsl_matrix_get_col(b, KS, i);
        gsl_linalg_QRPT_solve(temp3, tau, permutation, b, x);
        gsl_matrix_set_row(KP, i, x);
    }


    //gsl_blas_dtrsm(CblasRight,CblasLower,CblasNoTrans,CblasNonUnit,1.0,temp3,KP);

    // Scaling KI with b0
    gsl_matrix_scale(KP, b0);

    /*

    // Saving the current matrix of KP into KP_current for calculations
    gsl_matrix* KP_current = gsl_matrix_alloc(12, 3); // = get_Matrix(14);
    gsl_matrix_memcpy(KP_current, KP);

    // Saving the current subtrahend into sub for calculations
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,KP, D_op, 0.0, temp1);
    //gsl_matrix_mul_elements(KP_current, D_op);
    gsl_matrix* sub = gsl_matrix_alloc(12, 12);
    gsl_matrix_memcpy(sub, eye);
    gsl_matrix_sub(sub, temp1);

    // Calculating the Division of sub and KP_current
    //KP = KP_current;
    gsl_blas_dtrsm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1.0,sub,KP);
    gsl_linalg_QRPT_decomp(sub, tau, permutation, &signum, norm);
    for (j = 0; j < KP->size1; j++) {
        gsl_matrix_get_col(b, KP, i);
        gsl_linalg_QRPT_solve(sub, tau, permutation, b, x);
        gsl_matrix_set_row(KP, i, x);
    }
    //gsl_matrix_div_elements(KP, KP_current);

     */
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
void calculate_AG(void){

    // Storing the Matrix B into the new Matrix B_current for further calculations
    //gsl_matrix* B_current = get_Matrix(3);

    // Copying the matrix B_current in B_current2 for second calculation
    //gsl_matrix* B_current2 = BSystem;

    // Storing the matrix KI in new matrix KI_current
    //gsl_matrix* KI_current = get_Matrix(5);

    // Computing temp_calc1: -->  temp_calc1 = - BSystem * KI
    gsl_matrix* temp_calc1 = gsl_matrix_alloc(12, 3);
    gsl_matrix* temp_calc3 = gsl_matrix_alloc(12, 12);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, BSystem, KI, 0.0, temp_calc1);

    //gsl_matrix_scale(temp_calc1, -1);           // Multiply temp_calc1 with -1
    //gsl_matrix* A_current = get_Matrix(1);       // Storing the Matrix A in new Matrix A_current
    //gsl_matrix* KP_current = get_Matrix(14);     // Storing the Matrix KP in new Matrix KP_current
    //gsl_matrix* C_current = C_op;       // Storing the Matrix C_current in new Matrix C_current
    gsl_matrix* A_current = gsl_matrix_alloc(12, 12);
    gsl_matrix_memcpy(A_current, ASystem);

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
    gsl_matrix* temp_calc2 = gsl_matrix_alloc(12, 12);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, KP, C_op, 0.0, temp_calc2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, BSystem, temp_calc2, 0.0, temp_calc3);
    gsl_matrix_sub(A_current, temp_calc3);
    // TODO check it --> will probably not work

    // Copying the System matrix L (12 x 12) into top left corner of AG (15 x 15)
    for (size_t i = 0; i < 12; i++) {
        for (size_t j = 0; j < 12; j++) {
            gsl_matrix_set(Ag,i, j, gsl_matrix_get(A_current, i, j));
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
void calculate_EWG(void){
    //EW_G = Ag;
    gsl_matrix_memcpy(EW_G, Ag);
    gsl_eigen_nonsymm(EW_G, eigenvalue3, workspace);
}
// TODO check this methd, pretty same like calculate_EWI

/*
 * Comment
 */
void matrix_Calculator_EWG(void){
    // Counter for termination of while loop
    size_t i = 0;

    // Pointer on vector_view containing all real parts of complex vector eigenvalue3
    a = gsl_vector_complex_real(eigenvalue3);

    // Storing the structural element vector of vector_view a into gsl_vector real
    //gsl_vector* real = &a->vector;

    // Finding the biggest value of vector real
    double max = gsl_vector_max(&a.vector);

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
void tune_matrix_EWG(void){

    // Initialize matrix EW_I1 with ground matrix EW_I
    EW_G1 = EW_G;
    gsl_eigen_nonsymm(EW_G1,eigenvalue4,workspace);

    // counter for while loop
    size_t j = 0;
    size_t i = 1;

    // The biggest real number of complex eigenvalue vector of EW_G
    gsl_vector_view saving = gsl_vector_complex_real(eigenvalue3);
    double real1 = gsl_vector_max(&saving.vector);
    //for (size_t k = 0; k < 15; k++) {
      //  printf("%f + i%f %e \n",eigenvalue3->data[k*eigenvalue3->stride], gsl_vector_complex_get(eigenvalue3,k), gsl_vector_get(&saving.vector,k));
    //}

    // The biggest real number of complex eigenvalue vector of EW_G1 which represents the current state
    gsl_vector_view saving2 = gsl_vector_complex_real(eigenvalue4);
    double real2 = gsl_vector_max(&saving2.vector);

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
void tune_KP(void){
    double b_end = B0;
    EWG_output = EW_G1;
    calculate_KP(KS, b_end);

}
// TODO check

/*
 * This method returns the saved values A0 and B0
 */
double scalar(size_t n){
    switch (n){
        case 1:
            return A0;
        case 2:
            return B0;
        default:
            return 0;
    }
}

/*
 * This Function executes all methods written above to calculate the
 * PI state space controller
 */
void calculating_PI_Controller(size_t zaehler){
    getInputParameter();
    update_BSystem(zaehler);
    changing_engineSpeed(savingMatrix(30), zaehler);
    calculate_Cop();
    calculate_Dop();
    calculate_KS();
    calculate_KI(KS, scalar(1));
    //calculate_Ai();
    //calculate_EWI();
    //matrix_Calculator_EWI();
    //tune_matrix_EWI();
    //calculate_KP(get_Matrix(4), scalar(2));
    //calculate_AG();
    //calculate_EWG();
    //matrix_Calculator_EWG();
    //tune_matrix_EWG();
    //tune_KP();
    rechanging_engineSpeed();
}

/*
 * This method returns the vectors C_in and D_in
 * For Debugging purposes
 */
gsl_vector* returnCinDin(size_t n){
    if (n == 0) {
        return C_in;
    }
    if (n == 1) {
        return D_in;
    } else
        return 0;
}
// TODO check if it works correctly and how performance change

void update_BSystem(size_t zaehler){

    double T_n1 = gsl_matrix_get(savingMatrix(31), zaehler, 0);
    double T_n2 = gsl_matrix_get(savingMatrix(31), zaehler, 1);
    double T_n3 = gsl_matrix_get(savingMatrix(31), zaehler, 2);
    gsl_matrix_set(BSystem, 6, 3, (KONSTANTE / (T_n1*T_n1)));
    gsl_matrix_set(BSystem, 8, 4, (KONSTANTE / (T_n2*T_n2)));
    gsl_matrix_set(BSystem, 10, 5, (KONSTANTE / (T_n3*T_n3)));

}

/*
 * Undo the values of system matrix A and Ainv
 * to prior values after calculations
 */
void rechanging_engineSpeed(void){

    gsl_matrix_set(ASystem, 6, 6, -7.0356);           // Setting value in Matrix A_7_7
    gsl_matrix_set(ASystem, 6, 7, -183.0609);         // Setting value in Matrix A_7_8
    gsl_matrix_set(ASystem, 8, 8, -7.0356);           // Setting value in Matrix A_9_9
    gsl_matrix_set(ASystem, 8, 9, -183.0609);         // Setting value in Matrix A_9_10
    gsl_matrix_set(ASystem, 10, 10, -7.0356);         // Setting value in Matrix A_11_11
    gsl_matrix_set(ASystem, 10, 11, -183.0609);       // Setting value in Matrix A_11_12
    gsl_matrix_set(A_Inv, 7, 6, -0.0055);       // Setting value in Matrix A_inv 8_7
    gsl_matrix_set(A_Inv, 7, 7, -0.0384);       // Setting value in Matrix A_inv 8_8
    gsl_matrix_set(A_Inv, 9, 8, -0.0055);       // Setting value in Matrix A_inv 10_9
    gsl_matrix_set(A_Inv, 9, 9, -0.0384);       // Setting value in Matrix A_inv 10_10
    gsl_matrix_set(A_Inv, 11, 10, -0.0055);     // Setting value in Matrix A_inv 12_11
    gsl_matrix_set(A_Inv, 11, 11, -0.0384);     // Setting value in Matrix A_inv 12_12
}
