//
// Created by Farid Bonakdar on 04.08.16.
//

#include <gsl/gsl_matrix.h>
#include <memory.h>
#include "InitTest.h"
#include "horizontalModel.h"
#include <string.h>
#include <complex.h>


/*
 * This Class is generated to initialize all data for
 * Benchmarking and testing.
 * The methods should run once in the beginning of the simulation
 * to obtain all necessary vectors and matrices.
 */
#if DEBUGGER == 1
FILE* Cmat;
FILE* Dmat;
FILE* alphar;
FILE* alphax;
FILE* alphay;
FILE* axFile;
FILE* ayFile;
FILE* betaFile;
FILE* FzFile;
FILE* FxFile;
FILE* FyFile;
FILE* muFile;
FILE* muxFile;
FILE* muyFile;
FILE* psiPPFile;
FILE* vFile;
FILE* vrFile;
FILE* stepUFile;
FILE* stepXFile;
#endif
FILE* DrehzahlenFile;
FILE* statLenkFile;
FILE* dynLenkFile;
FILE* xgFile;
FILE* ugFile;
FILE* accFile;
FILE* accAltFile;
#if DEBUGGER == 1
FILE* axtermeFile;
FILE* aytermeFile;
FILE* acc2File;
FILE* xg_diff_loopFile;
FILE* xg_alt_diffFile;
FILE* ug_diff_loopFile;
FILE* ug_alt_diffFile;
FILE* acc_diff_loopFile;
#endif
//FILE* EWgFile;
//FILE* EWiFile;
//FILE* EWgAnfangFile;
//FILE* EWgwhile1File;
//FILE* EWiAnfangFile;
//FILE* EWiwhile1File;
FILE* KiFile;
FILE* KiAnfangFile;
FILE* KiEndeFile;
FILE* KpFile;
FILE* KpAnfangFile;
FILE* KpEndeFile;
FILE* KptransFile;
FILE* KsOut1File;
FILE* AgAnfangFile;
FILE* AiOutFile;
FILE* n_updnFile;
FILE* T_nFile;
FILE* prod1File;
FILE* prod2File;

gsl_vector* n1;
gsl_vector* n2;
gsl_vector* n3;
gsl_vector* delta1;
gsl_vector* delta2;
gsl_vector* delta3;
gsl_vector* deltaDyn1;
gsl_vector* deltaDyn2;
gsl_vector* deltaDyn3;
gsl_vector* velocity_x;
gsl_vector* velocity_y;
gsl_vector* yawrate;
gsl_matrix* cMatrix;
gsl_matrix* dMatrix;
gsl_vector* a_x;
gsl_vector* a_y;
gsl_matrix* alphaRMatrix;
gsl_matrix* alphaXMatrix;
gsl_matrix* alphaYMatrix;
gsl_matrix* betaMatrix;
gsl_matrix* FxMatrix;
gsl_matrix* FyMatrix;
gsl_matrix* FzMatrix;
gsl_matrix* muMatrix;
gsl_matrix* muXMatrix;
gsl_matrix* muYMatrix;
gsl_vector* psiPPVec;
gsl_matrix* stepXMatrix;
gsl_matrix* stepUMatrix;
gsl_matrix* vMatrix;
gsl_matrix* vrMatrix;
gsl_matrix* DrehzahlenMat;
gsl_matrix* statLenkMatrix;
gsl_matrix* dynLenkMatrix;
gsl_matrix* xgMatrix;
gsl_matrix* ugMatrix;
gsl_matrix* accMatrix;
gsl_matrix* accAltMatrix;
gsl_matrix* axtermeMatrix;
gsl_matrix* aytermeMatrix;
gsl_matrix* acc2Matrix;
gsl_matrix* xg_diff_loopMatrix;
gsl_matrix* xg_alt_diffMatrix;
gsl_matrix* ug_diff_loopMatrix;
gsl_matrix* ug_alt_diffMatrix;
gsl_matrix* acc_diff_loopMatrix;
float gesamtzeit = 0.0;
//gsl_matrix_complex* EWgMatrix;
//gsl_matrix_complex* EWiMatrix;
//gsl_matrix_complex* EWgAnfangMatrix;
//gsl_matrix_complex* EWgwhile1Matrix;
//gsl_matrix_complex* EWiAnfangMatrix;
//gsl_matrix_complex* EWiwhile1Matrix;
gsl_matrix* KiMatrix;
gsl_matrix* KiAnfangMatrix;
gsl_matrix* KiEndeMatrix;
gsl_matrix* KpMatrix;
gsl_matrix* KpAnfangMatrix;
gsl_matrix* KpEndeMatrix;
gsl_matrix* KptransMatrix;
gsl_matrix* KsOut1Matrix;
gsl_matrix* AgAnfangMatrix;
gsl_matrix* AiOutMatrix;
gsl_vector* tempCol;
gsl_matrix* n_updnMatrix;
gsl_matrix* T_nMatrix;
gsl_matrix* prod1Matrix;
gsl_matrix* prod2Matrix;



/*
 * This method initialize all vectors
 * The size of the vector must be the same as the file size
 */
void start_initializing(size_t choice) {
    if(choice == 0){
        n1 = gsl_vector_alloc(61001);
        n2 = gsl_vector_alloc(61001);
        n3 = gsl_vector_alloc(61001);
        delta1 = gsl_vector_alloc(61001);
        delta2 = gsl_vector_alloc(61001);
        delta3 = gsl_vector_alloc(61001);
        deltaDyn1 = gsl_vector_alloc(61001);
        deltaDyn2 = gsl_vector_alloc(61001);
        deltaDyn3 = gsl_vector_alloc(61001);
        velocity_x = gsl_vector_alloc(61001);
        velocity_y = gsl_vector_alloc(61001);
        yawrate = gsl_vector_alloc(61001);
        cMatrix = gsl_matrix_alloc(61001, 18);
        dMatrix = gsl_matrix_alloc(61001, 18);
        a_x = gsl_vector_alloc(61001);
        a_y = gsl_vector_alloc(61001);
        alphaRMatrix = gsl_matrix_alloc(61001, 3);
        alphaXMatrix = gsl_matrix_alloc(61001, 3);
        alphaYMatrix = gsl_matrix_alloc(61001, 3);
        betaMatrix = gsl_matrix_alloc(61001, 3);
        FxMatrix = gsl_matrix_alloc(61001, 3);
        FyMatrix = gsl_matrix_alloc(61001, 3);
        FzMatrix = gsl_matrix_alloc(61001, 3);
        muMatrix = gsl_matrix_alloc(61001, 3);
        muXMatrix = gsl_matrix_alloc(61001, 3);
        muYMatrix = gsl_matrix_alloc(61001, 3);
        psiPPVec = gsl_vector_alloc(61001);
        stepXMatrix = gsl_matrix_alloc(61001, 3);
        stepUMatrix = gsl_matrix_alloc(61001, 9);
        vMatrix = gsl_matrix_alloc(61001, 3);
        vrMatrix = gsl_matrix_alloc(61001, 3);
        DrehzahlenMat = gsl_matrix_alloc(61001, 3);
        statLenkMatrix = gsl_matrix_alloc(61001, 3);
        dynLenkMatrix = gsl_matrix_alloc(61001, 3);
        xgMatrix = gsl_matrix_alloc(61001, 3);
        ugMatrix = gsl_matrix_alloc(61001, 9);
        accMatrix = gsl_matrix_alloc(61001, 3);
        accAltMatrix = gsl_matrix_alloc(61001, 3);
#if DEBUGGER == 1
        axtermeMatrix = gsl_matrix_alloc(61001, 7);
        aytermeMatrix = gsl_matrix_alloc(61001, 4);
        acc2Matrix = gsl_matrix_alloc(61001, 3);
        acc_diff_loopMatrix = gsl_matrix_alloc(61001, 3);
        xg_diff_loopMatrix = gsl_matrix_alloc(61001, 3);
        xg_alt_diffMatrix = gsl_matrix_alloc(61001, 3);
        ug_diff_loopMatrix = gsl_matrix_alloc(61001, 9);
        ug_alt_diffMatrix = gsl_matrix_alloc(61001, 9);
#endif

        // Data for linear.c
        //EWgMatrix = gsl_matrix_complex_alloc(61001, 15);
        //EWiMatrix = gsl_matrix_complex_alloc(61001, 15);
        //EWgAnfangMatrix = gsl_matrix_complex_alloc(61001, 15);
        //EWgwhile1Matrix = gsl_matrix_complex_alloc(61001, 15);
        //EWiAnfangMatrix = gsl_matrix_complex_alloc(61001, 15);
        //EWiwhile1Matrix = gsl_matrix_complex_alloc(61001, 15);
        KiMatrix = gsl_matrix_alloc(36, 61001);
        KiAnfangMatrix = gsl_matrix_alloc(36, 61001);
        KiEndeMatrix = gsl_matrix_alloc(36, 61001);
        //KpMatrix = gsl_matrix_alloc(36, 61001);
        //KpAnfangMatrix = gsl_matrix_alloc(36, 61001);
        //KpEndeMatrix = gsl_matrix_alloc(36, 61001);
        //KptransMatrix = gsl_matrix_alloc(36, 61001);
        KsOut1Matrix = gsl_matrix_alloc(36, 61001);
        //AgAnfangMatrix = gsl_matrix_alloc(225, 61001);
        //AiOutMatrix = gsl_matrix_alloc(225, 61001);
        tempCol = gsl_vector_alloc(9);
        n_updnMatrix = gsl_matrix_alloc(61001, 3);
        T_nMatrix = gsl_matrix_alloc(61001, 3);
        prod1Matrix = gsl_matrix_alloc(9, 61001);
        prod2Matrix = gsl_matrix_alloc(9, 61001);
    }
}

/*
 * Method starts to scan the values of each file into its vector
 * and saves it
 */
void start_reading(void) {
    gsl_matrix_get_col(n1, DrehzahlenMat, 0);
    gsl_matrix_get_col(n2, DrehzahlenMat, 1);
    gsl_matrix_get_col(n3, DrehzahlenMat, 2);
    gsl_matrix_get_col(delta1, statLenkMatrix, 0);
    gsl_matrix_get_col(delta2, statLenkMatrix, 1);
    gsl_matrix_get_col(delta3, statLenkMatrix, 2);
    gsl_matrix_get_col(deltaDyn1, dynLenkMatrix, 0);
    gsl_matrix_get_col(deltaDyn2, dynLenkMatrix, 1);
    gsl_matrix_get_col(deltaDyn3, dynLenkMatrix, 2);

#if DEBUGGER == 1
    gsl_matrix_fscanf(Cmat, cMatrix);
    gsl_matrix_fscanf(Dmat, dMatrix);
    gsl_matrix_fscanf(alphar, alphaRMatrix);
    gsl_matrix_fscanf(alphax, alphaXMatrix);
    gsl_matrix_fscanf(alphay, alphaYMatrix);
    gsl_vector_fscanf(axFile, a_x);
    gsl_vector_fscanf(ayFile, a_y);
    gsl_matrix_fscanf(betaFile, betaMatrix);
    gsl_matrix_fscanf(FxFile, FxMatrix);
    gsl_matrix_fscanf(FyFile, FyMatrix);
    gsl_matrix_fscanf(FzFile, FzMatrix);
    gsl_matrix_fscanf(muFile, muMatrix);
    gsl_matrix_fscanf(muyFile, muYMatrix);
    gsl_matrix_fscanf(muxFile, muXMatrix);
    gsl_vector_fscanf(psiPPFile, psiPPVec);
    gsl_matrix_fscanf(vFile, vMatrix);
    gsl_matrix_fscanf(vrFile, vrMatrix);
    gsl_matrix_fscanf(stepXFile, stepXMatrix);
    gsl_matrix_fscanf(stepUFile, stepUMatrix);
#endif
    gsl_matrix_fscanf(DrehzahlenFile, DrehzahlenMat);
    gsl_matrix_fscanf(statLenkFile, statLenkMatrix);
    gsl_matrix_fscanf(dynLenkFile, dynLenkMatrix);
    gsl_matrix_fscanf(xgFile, xgMatrix);
    gsl_matrix_fscanf(ugFile, ugMatrix);
    gsl_matrix_fscanf(accFile, accMatrix);
    gsl_matrix_fscanf(accAltFile, accAltMatrix);
#if DEBUGGER == 1
    gsl_matrix_fscanf(axtermeFile, axtermeMatrix);
    gsl_matrix_fscanf(aytermeFile, aytermeMatrix);
    gsl_matrix_fscanf(acc2File, acc2Matrix);
    gsl_matrix_fscanf(acc_diff_loopFile, acc_diff_loopMatrix);
    gsl_matrix_fscanf(xg_diff_loopFile, xg_diff_loopMatrix);
    gsl_matrix_fscanf(xg_alt_diffFile, xg_alt_diffMatrix);
    gsl_matrix_fscanf(ug_diff_loopFile, ug_diff_loopMatrix);
    gsl_matrix_fscanf(ug_alt_diffFile, ug_alt_diffMatrix);
#endif

    // Matrices for linear.c
    //gsl_matrix_complex_fscanf(EWgFile, EWgMatrix);
    //gsl_matrix_complex_fscanf(EWgAnfangFile, EWgAnfangMatrix);
    //gsl_matrix_complex_fscanf(EWgwhile1File, EWgwhile1Matrix);
    //gsl_matrix_complex_fscanf(EWiFile, EWiMatrix);
    //gsl_matrix_complex_fscanf(EWiAnfangFile, EWiAnfangMatrix);
    //gsl_matrix_complex_fscanf(EWiwhile1File, EWiwhile1Matrix);
    //gsl_matrix_fscanf(EWiwhile1File, EWiwhile1Matrix);
    gsl_matrix_fscanf(KiFile, KiMatrix);
    gsl_matrix_fscanf(KiAnfangFile, KiAnfangMatrix);
    gsl_matrix_fscanf(KiEndeFile, KiEndeMatrix);
    //gsl_matrix_fscanf(KpFile, KpMatrix);
    //gsl_matrix_fscanf(KpAnfangFile, KpAnfangMatrix);
    //gsl_matrix_fscanf(KpEndeFile, KpEndeMatrix);
    //gsl_matrix_fscanf(KptransFile, KptransMatrix);
    gsl_matrix_fscanf(KsOut1File, KsOut1Matrix);
    //gsl_matrix_fscanf(AgAnfangFile, AgAnfangMatrix);
    //gsl_matrix_fscanf(AiOutFile, AiOutMatrix);
    gsl_matrix_fscanf(n_updnFile, n_updnMatrix);
    gsl_matrix_fscanf(T_nFile, T_nMatrix);
    gsl_matrix_fscanf(prod1File, prod1Matrix);
    gsl_matrix_fscanf(prod2File, prod2Matrix);
}

/*
 * This method gets an input number and returns a vector
 * which will be stored for later calculations
 */
gsl_vector* saving(size_t n) {
    switch (n){
        case 1:
            return n1;
        case 2:
            return n2;
        case 3:
            return n3;
        case 4:
            return delta1;
        case 5:
            return delta2;
        case 6:
            return delta3;
        case 7:
            return deltaDyn1;
        case 8:
            return deltaDyn2;
        case 9:
            return deltaDyn3;
        case 10:
            return velocity_x;
        case 11:
            return velocity_y;
        case 12:
            return yawrate;
        case 13:
            return a_x;
        case 14:
            return a_y;
        case 15:
            return psiPPVec;
        default:
            return 0;
    }
}

/*
 * This Method gets an input number and returns a matrix
 * which will be stored for later calculations
 */
gsl_matrix* savingMatrix(size_t n){
    switch (n){
        case 1:
            return cMatrix;
        case 2:
            return dMatrix;
        case 3:
            return alphaRMatrix;
        case 4:
            return alphaXMatrix;
        case 5:
            return alphaYMatrix;
        case 6:
            return betaMatrix;
        case 7:
            return FxMatrix;
        case 8:
            return FyMatrix;
        case 9:
            return FzMatrix;
        case 10:
            return muMatrix;
        case 11:
            return muXMatrix;
        case 12:
            return muYMatrix;
        case 13:
            return stepXMatrix;
        case 14:
            return stepUMatrix;
        case 15:
            return vMatrix;
        case 16:
            return vrMatrix;
        case 17:
            return DrehzahlenMat;
        case 18:
            return xgMatrix;
        case 19:
            return ugMatrix;
        case 20:
            return accMatrix;
        case 21:
            return accAltMatrix;
        case 22:
            return axtermeMatrix;
        case 23:
            return aytermeMatrix;
        case 24:
            return acc2Matrix;
        case 25:
            return acc_diff_loopMatrix;
        case 26:
            return xg_diff_loopMatrix;
        case 27:
            return xg_alt_diffMatrix;
        case 28:
            return ug_diff_loopMatrix;
        case 29:
            return ug_alt_diffMatrix;
        case 30:
            return n_updnMatrix;
        case 31:
            return T_nMatrix;
        default:
            return 0;
    }
}

/*
 *
 */
gsl_matrix* bigMatrices(size_t n){
    switch (n){
        case 1:
            return KsOut1Matrix;
        case 2:
            return KiMatrix;
        case 3:
            return KiAnfangMatrix;
        case 4:
            return prod1Matrix;
        case 5:
            return prod2Matrix;
        default:
            return 0;
    }

}

/*
 * This Method gets an input counter, a matrix pointer to read from and a pointer to
 * save the information.
 * The output is the transformation of a colon vector from matrixReader to a shaped
 * matrixSaver
 *      __________________________________
 *     | x11  x12  x13  x14  x15 ...  x1m |
 *     | x21  x22  x23  x24  x25 ...  x2m |
 *     | x31  x32  ...                    |
 * --> | x41  ...                      .  |
 *     | x51                           .  |     <---
 *     |  .        ...                    |     <---    Big Matrix with all test Data
 *     |  .                  ...          |     <---
 *     |  .                               |
 *     | xn1  ...                     xnm |
 *     |  |    |    |                  |  |
 *     |__|____|____|__________________|__|
 *
 *        1    2    3   ...            m        <---    Colon Vectors to read
 *
 *
 *        1:
 *      ___________________________________
 *     | x11  x41  x71  ...          xn-21 |
 *     | x21  x51  x81  ...          xn-11 |    <---    Matrix build with the colon vector
 *     | x31  x61  x91  ...           xn1  |    <---    with the shape of 3 x 12 for example
 *     |___________________________________|
 *
 */
void create_data(size_t zaehler, gsl_matrix* matrixReader, gsl_matrix* matrixSaver) {
    // Get the current colon from the big matrixReader and save it to temporary vector
    gsl_matrix_get_col(tempCol, matrixReader, zaehler);

    // Save the desired shape matrixSaver for further fill
    size_t rows = matrixSaver->size1;
    size_t cols = matrixSaver->size2;
    size_t total = matrixReader->size1;

    // Declare counter integers
    size_t i;
    size_t j;
    size_t k = 0;

    // For i is smaller then current cols of matrixSaver
    for (i = 0; i < cols ; i++) {
        // For j is smaller the current rows of matrixSaver
        for (j = 0; j < rows; j++) {
            // If k have not reached last element of colon vector of matrixReader do
            if (k!=total) {
                // Get value of temporary vector
                double zahl;
                zahl = gsl_vector_get(tempCol, k);
                // Save value of temporary vector to matrixSaver at position i x j
                gsl_matrix_set(matrixSaver, j, i, zahl);

                // Increment k for position
                k++;
            } else
                break;
        }
    }


}

/*
 * This Method is written for debugging process.
 * It gets a matrix or a vector an show all values of
 * the structure in the terminal
 */
void printer(gsl_matrix* matrix, gsl_vector* vector){
    if (matrix != NULL && vector == NULL) {
        size_t rows = matrix->size1;
        size_t col = matrix->size2;
        size_t i;
        size_t j;
        printf("_________________");
        for (i = 0; i < rows; i++) {
            printf(" \n");
            for (j = 0; j < col; j++) {
                printf("%.16f ", gsl_matrix_get(matrix, i, j));
            }
        }
        printf("\n_________________");
    }
    if (matrix == NULL && vector != NULL){
        size_t rows = vector->size;
        size_t i;
        printf("[ ");
        for (i = 0; i < rows; i++) {
            printf("%.6f ",gsl_vector_get(vector, i));
        }
        printf(" ]\n");
    }
}

/*
 * This method opens all the txt files in directory "Raw data"
 * and saves them into its FILE Container
 *
 *
 */
void open_files(size_t choice) {
    /*
     * Change your current path for testing
     * Beware of Operating System and use of slash
     * and type in the full path of your data
     * OSX/Linux --> /
     * Windows   --> \
     */
    if(choice == 0){
#if DEBUGGER == 1
        Cmat = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/CMatrix.txt", "rw");
        Dmat = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/DMatrix.txt", "rw");
        alphar = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/alpha_r.txt", "rw");
        alphax = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/alpha_x.txt", "rw");
        alphay = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/alpha_y.txt", "rw");
        axFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/ax.txt", "rw");
        ayFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/ay.txt", "rw");
        betaFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/beta.txt", "rw");
        FzFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/Fz.txt", "rw");
        FxFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/Fx.txt", "rw");
        FyFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/Fy.txt", "rw");
        muFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/mu.txt", "rw");
        muxFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/mux.txt", "rw");
        muyFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/muy.txt", "rw");
        psiPPFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/psi_pp.txt", "rw");
        vFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/Geschwindigkeit.txt", "rw");
        vrFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/vr.txt", "rw");
        stepUFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/step_u.txt", "rw");
        stepXFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/step_x.txt", "rw");
#endif
        DrehzahlenFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/Drehzahl.txt", "rw");
        statLenkFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/statischer_Lenkwinkel.txt", "rw");
        dynLenkFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/dynLenk.txt", "rw");
        xgFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/xg.txt", "rw");
        ugFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/ug.txt", "rw");
        accFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/acc.txt", "rw");
        accAltFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/acc_alt.txt", "rw");

#if DEBUGGER == 1
        axtermeFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/ax_terme.txt", "rw");
        aytermeFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/ay_terme.txt", "rw");
        acc2File = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/acc2.txt", "rw");
        xg_diff_loopFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/xg_diff_loop.txt", "rw");
        xg_alt_diffFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/xg_alt_diff.txt", "rw");
        ug_diff_loopFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/ug_diff_loop.txt", "rw");
        ug_alt_diffFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/ug_alt_diff.txt", "rw");
        acc_diff_loopFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/acc_diff_loop.txt", "rw");
#endif

        //EWgFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/EWgNeu.txt", "rw");
        //EWiFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/EWiNeu.txt", "rw");
        //EWgAnfangFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/EWganfang.txt", "rw");
        //EWgwhile1File = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/EWgwhile1.txt", "rw");
        //EWiAnfangFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/EWianfang.txt", "rw");
        //EWiwhile1File = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/EWiwhile1.txt", "rw");
        KiFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/KiNeu.txt", "rw");
        KiAnfangFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/KianfangNeu.txt", "rw");
        KiEndeFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/KiendeNeu.txt", "rw");
        //KpFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/KpNeu.txt", "rw");
        //KpAnfangFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/KpanfangNeu.txt", "rw");
        //KpEndeFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/KpendeNeu.txt", "rw");
        //KptransFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/KptransNeu.txt", "rw");
        KsOut1File = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/Ksout1Neu.txt", "rw");
        //AgAnfangFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/AganfangNeu.txt", "rw");
        //AiOutFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/Aiout1Neu.txt", "rw");
        n_updnFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/n_updn_neu.txt", "rw");
        T_nFile = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/T_n.txt", "rw");
        prod1File = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/prod1KI.txt", "rw");
        prod2File = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Datenanalyse/INVNEU.txt", "rw");
    }
}

/*
 * This method gets an array of cycle times an prints the
 * desired output.
 * --> "MAX" = longest calculation time
 * --> "MIN" = shortest calculation time
 * --> "TOTAL" = whole execution time
 * e.g
 * .... calculateCycleTime(timings, "MAX");
 */
void calculateCycleTime (float timings[], char *auswahl){
    // Initializing the output times with random values
    float maxTime = timings[0];

    // Minimal time must be big enough to get a sensible output
    float minTime = 200.0;

    // Starting with 0
    //float gesamtzeit = 0.0;

    // Declaration of counters used in for-loop
    size_t i, j, k;
    size_t count = sizeof(&timings)/ sizeof(float);

    /*
     * If auswahl is "MAX" iterate over all array entries and
     * save the bigger time in maxTime
     * else continue
     * At the end print the value
     */
    if(compare_strings("MAX", auswahl) == 0) {
        for (i = 0; i< count; i++) {
            if (maxTime <= timings[i]){
                maxTime = timings[i];
            } else continue;
        }
        printf(" Maximalzeit betraegt %g us\n", 1000000*maxTime);
    }

    /*
     * If auswahl is "MIN" iterate over all array entries and
     * save the smallest time in minTime
     * else continue
     * At the end print the value
     */
    if(compare_strings("MIN", auswahl) == 0) {
        for (j = 0; j< count; j++) {
            if (minTime >= timings[j]){
                minTime = timings[j];
            } else continue;
        }
        printf(" Minimalzeit betraegt %g us\n", 1000000*minTime);
    }

    /*
     * If auswahl is "Total" iterate over all array entries and
     * increment the value to gesamtzeit
     * At the end print the value
     */
    if(compare_strings("Total", auswahl) == 0) {
        for (k = 0; k< count; k++) {
                gesamtzeit += timings[k];
        }
        printf(" Gesamtzeit betraegt %g us\n", 1000000*gesamtzeit);
    }

}

/*
 * This Method gets two time clocks and calculates the execution time
 * in respect to the arch specific clocks per second and return it
 * as float
 */
float saveTiming(clock_t Anfang, clock_t Ende){
    float time = (float) (Ende - Anfang) / CLOCKS_PER_SEC;
    return time;
}

/*
 * This method compares the two input strings with each other
 * and return 0
 */
int compare_strings(char a[], char b[]) {

    size_t counter;
    counter = 0;

    // This while statement compares the lengths of the two strings
    while (a[counter] == b[counter]) {
        if (a[counter] == '\0' || b[counter] == '\0')
            break;
        counter++;
    }

    // If content of two strings is the same then return value is 0 else -1
    if (a[counter] == '\0' && b[counter] == '\0')
        return 0;
    else
        return -1;
}

/*
 * This method prints out how much time the calculation takes
 */
void print_Timings(float timings[],size_t size){
    size_t j;
    for (j = 0;j<size;j++){
        printf(" Zeit : %g\n",timings[j]);
    }

}

/*
 * TODO comment this function
 */
void drehzahlTester(void){
    gsl_matrix_get_col(n1, DrehzahlenMat, 0);
    gsl_matrix_get_col(n2, DrehzahlenMat, 1);
    gsl_matrix_get_col(n3, DrehzahlenMat, 2);
    gsl_matrix_get_col(delta1, statLenkMatrix, 0);
    gsl_matrix_get_col(delta2, statLenkMatrix, 1);
    gsl_matrix_get_col(delta3, statLenkMatrix, 2);
    gsl_matrix_get_col(deltaDyn1, dynLenkMatrix, 0);
    gsl_matrix_get_col(deltaDyn2, dynLenkMatrix, 1);
    gsl_matrix_get_col(deltaDyn3, dynLenkMatrix, 2);
}
