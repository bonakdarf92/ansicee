//
// Created by Farid Bonakdar on 04.08.16.
//

#include "InitTest.h"
#include "horizontalModel.h"

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
void start_initializing() {
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
void start_reading() {
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
    }
}

void create_data() {
    // TODO Implement core
}
void printer(){
}

/*
 * This method opens all the txt files in directory "Raw data"
 * and saves them into its FILE Container
 *
 *
 */
void open_files() {
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
