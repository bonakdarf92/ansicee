//
// Created by Farid Bonakdar on 04.08.16.
//

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <memory.h>
#include "InitTest.h"
#include "horizontalModel.h"
#include <string.h>


/*
 * This Class is generated to initialize all data for
 * Benchmarking and testing.
 * The methods should run once in the beginning of the simulation
 * to obtain all necessary vectors and matrices.
 */


FILE* n_1;
FILE* n_2;
FILE* n_3;
FILE* delta_1;
FILE* delta_2;
FILE* delta_3;
FILE* delta_dyn_1;
FILE* delta_dyn_2;
FILE* delta_dyn_3;
FILE* geschwindigkeit_x;
FILE* geschwindigkeit_y;
FILE* gierrate;
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



/*
 * This method initialize all vectors
 * The size of the vector must be the same as the file size
 */
void start_initializing(void) {
    n1 = gsl_vector_alloc(6001);
    n2 = gsl_vector_alloc(6001);
    n3 = gsl_vector_alloc(6001);
    delta1 = gsl_vector_alloc(6001);
    delta2 = gsl_vector_alloc(6001);
    delta3 = gsl_vector_alloc(6001);
    deltaDyn1 = gsl_vector_alloc(6001);
    deltaDyn2 = gsl_vector_alloc(6001);
    deltaDyn3 = gsl_vector_alloc(6001);
    velocity_x = gsl_vector_alloc(6001);
    velocity_y = gsl_vector_alloc(6001);
    yawrate = gsl_vector_alloc(6001);

}

/*
 * Method starts to scan the values of each file into its vector
 * and saves it
 */
void start_reading(void) {
    gsl_vector_fscanf(n_1, n1);
    gsl_vector_fscanf(n_2, n2);
    gsl_vector_fscanf(n_3, n3);
    gsl_vector_fscanf(delta_1, delta1);
    gsl_vector_fscanf(delta_2, delta2);
    gsl_vector_fscanf(delta_3, delta3);
    gsl_vector_fscanf(delta_dyn_1, deltaDyn1);
    gsl_vector_fscanf(delta_dyn_2, deltaDyn2);
    gsl_vector_fscanf(delta_dyn_3, deltaDyn3);
    gsl_vector_fscanf(geschwindigkeit_x, velocity_x);
    gsl_vector_fscanf(geschwindigkeit_y, velocity_y);
    gsl_vector_fscanf(gierrate, yawrate);

}

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
        default:
            return 0;
    }
}

void create_data(void) {
    // TODO Implement core
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
                printf("%f ", gsl_matrix_get(matrix, i, j));
            }
        }
        printf("\n_________________");
    }
    if (matrix == NULL && vector != NULL){
        size_t rows = vector->size;
        size_t i;
        printf("[ ");
        for (i = 0; i < rows; i++) {
            printf("%f ",gsl_vector_get(vector, i));
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
void open_files(void) {
    /*
     * Change your current path for testing
     * Beware of Operating System and use of slash
     * and type in the full path of your data
     * OSX/Linux --> /
     * Windows   --> \
     */
    n_1 = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Raw Data/n_1.txt", "rw");
    n_2 = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Raw data/n_2.txt", "rw");
    n_3 = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Raw data/n_3.txt", "rw");
    delta_1 = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Raw data/delta_1.txt", "rw");
    delta_2 = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Raw data/delta_2.txt", "rw");
    delta_3 = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Raw data/delta_3.txt", "rw");
    delta_dyn_1 = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Raw data/delta_dyn_1.txt", "rw");
    delta_dyn_2 = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Raw data/delta_dyn_2.txt", "rw");
    delta_dyn_3 = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Raw data/delta_dyn_3.txt", "rw");
    geschwindigkeit_x = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Raw data/geschwindigkeit_x.txt", "rw");
    geschwindigkeit_y = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Raw data/geschwindigkeit_y.txt", "rw");
    gierrate = fopen("/Users/faridbonakdar/ClionProjects/ansicee/Benchmark/Raw data/gierrate.txt", "rw");
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
    float gesamtzeit = 0.0;

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
