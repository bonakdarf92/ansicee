//
// Created by Farid Bonakdar on 30.06.16.
//

#include "calcMatrix_A_B.h"

/*
 * Hier werden die Berechnungsvorschriften fuer die Matrixoperationen geschrieben.
 * Auszufuehren sind mindestens die unten implementierten Funktionen.
 * Bei Bedarf evtl. eine extra Bibliothek inkludieren um performanz zu steigern.
 *
 */

void invertMatrixA(gsl_matrix* matrix, gsl_matrix* output) {
    size_t rows1 = matrix->size1;                       // Number of rows from Matrix matrix
    size_t cols1 = matrix->size2;                       // Number of cols from Matrix matrix
    size_t i, j;

    gsl_matrix* transMatrix;
    transMatrix = gsl_matrix_alloc(rows1, cols1);       // Allocating the transpose of Matrix matrix
    gsl_matrix_transpose_memcpy(transMatrix, matrix);   // Transpose Matrix matrix into transMatrix

    // Allocating Matrices and vector for singular value decomposition
    gsl_matrix* U = gsl_matrix_alloc(cols1, rows1);
    gsl_matrix* V = gsl_matrix_alloc(rows1, cols1);
    gsl_vector* S = gsl_vector_alloc(rows1);
    gsl_vector* work = gsl_vector_alloc(rows1);

    // SVD
    gsl_linalg_SV_decomp(transMatrix, V, S, work);
    gsl_vector_free(work);
    gsl_matrix_memcpy(U, transMatrix);

    // Inversion of Matrix S
    gsl_matrix* SP = gsl_matrix_alloc(rows1, rows1);
    gsl_matrix_set_zero(SP);
    for (i = 0; i < rows1; ++i) {
        gsl_matrix_set(SP, i, i, gsl_vector_get(S, i));
    }
    gsl_permutation* permutation = gsl_permutation_alloc(rows1);
    int sgn;
    gsl_linalg_LU_decomp(SP, permutation, &sgn);

    gsl_matrix* SI = gsl_matrix_calloc(rows1, rows1);

    for (i = 0; i < rows1; ++i) {
        if (gsl_vector_get(S, i) > 0.000001)
            gsl_matrix_set(SI, i, i, 1.0 / gsl_vector_get(S, i));
    }

    gsl_matrix* VT = gsl_matrix_alloc(rows1, rows1);
    gsl_matrix_transpose_memcpy(VT, V);


    // PseudoInverse
    // -------------------------------------
    // +++++++++++++++++++++++++++++++++++++
    // -------------------------------------

    gsl_matrix* SIpVT = gsl_matrix_alloc(rows1, rows1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, SI, VT, 0.0, SIpVT);
    gsl_matrix* pInv = gsl_matrix_alloc(cols1, rows1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, SIpVT, 0.0, pInv);

    // End Pseudoinverse
    // _____________________________________
    gsl_matrix_memcpy(output, pInv);
    /*
    gsl_matrix_free(VT);
    gsl_matrix_free(SI);
    gsl_matrix_free(SIpVT);
    gsl_matrix_free(transMatrix);
    gsl_matrix_free(U);
    gsl_vector_free(S);
    gsl_matrix_free(V);
*/
}

void invertMatrixB(double (*Matrix)[12]) {

}

double eigenWerte() {
    return 0;
}

double determinante(gsl_matrix* matrix) {
    size_t rows = matrix->size1;
    size_t col = matrix->size2;
    if (rows != col){
        printf("");
        return 0;
    }
    return 0;
}













