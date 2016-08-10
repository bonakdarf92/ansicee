//
// Created by Farid Bonakdar on 30.06.16.
//

#include "horizontalModel.h"
#include <math.h>
#include "InitTest.h"


/*
 * Declaration of System matrices and inner System vector for further use in methods.
 * Stored here to get access to each scalar, vector and matrix.
 */

    gsl_vector* xg;       // Vector xg
    gsl_vector* ug;       // Vector ug
    gsl_vector* xg_alt;   // Vector xg_alt
    gsl_vector* ug_alt;   // Vector ug_Alt
    gsl_vector* acc_alt;  // Vector acc_alt
    gsl_vector* delta_x;  // Vector delta_x
    gsl_vector* delta_u;  // Vector delta_u
    gsl_vector* alpha_r;  // Vector alpha_r
    gsl_vector* alpha_x;   // Vector alpha_x
    gsl_vector* alpha_y;   // Vector alpha_y
    gsl_vector* C;        // Vector C
    gsl_vector* D;        // Vector D
    double v_1;           // Velocity No.1
    double v_2;           // Velocity No.2
    double v_3;           // Velocity No.3
    double beta_1;        // Angle beta_1
    double beta_2;        // Angle beta_2
    double beta_3;        // Angle beta_3
    gsl_vector* beta;     // Vector beta
    gsl_vector* v;        // Vector v
    gsl_vector* vr;       // Vector vr
    gsl_vector* sr;       // Vector sr
    gsl_vector* mu;       // Vector mu
    gsl_vector* mux;      // Vector mu
    gsl_vector* muy;      // Vector mu
    gsl_vector* Fz;       // Force Fz
    gsl_vector* Fx;       // Force Fx
    gsl_vector* Fy;       // Force Fy
    gsl_vector* acc;      // Vector acc
    gsl_matrix* test_ug;  // Matrix for testing and simulation ug
    gsl_matrix* test_xg;  // Matrix for testing and simulation xg


/*
 * Initialization of System matrix C and D.
 * Allocation of inner System vector xg, ug, xg_alt, ug_alt and acc_alt
 * The size of vector and matrix taken from Jan's Matlab Code (Horizontal model)
 */
void initializeVector(){

    xg = gsl_vector_alloc(3);       // Vector xg
    ug = gsl_vector_alloc(9);       // Vector ug
    xg_alt = gsl_vector_alloc(3);   // Vector xg_alt
    ug_alt = gsl_vector_alloc(9);   // Vector ug_Alt
    acc_alt = gsl_vector_alloc(3);  // Vector acc_alt
    delta_x = gsl_vector_alloc(3);  // Vector delta_x
    delta_u = gsl_vector_alloc(9);  // Vector delta_u
    alpha_r = gsl_vector_alloc(3);  // Vector alpha_r
    alpha_x = gsl_vector_alloc(3);   // Vector alpha_x
    alpha_y = gsl_vector_alloc(3);   // Vector alpha_y
    C = gsl_vector_alloc(18);    // Matrix/Vector C
    D = gsl_vector_alloc(18);    // Matrix/Vector D
    beta = gsl_vector_alloc(3);     // Vector beta
    v = gsl_vector_alloc(3);        // Vector v
    vr = gsl_vector_alloc(3);       // Vector vr
    sr = gsl_vector_alloc(3);       // Vector sr
    mu = gsl_vector_alloc(3);       // Vector mu
    mux = gsl_vector_alloc(3);      // Vector mux
    muy = gsl_vector_alloc(3);      // Vector muy
    Fz = gsl_vector_alloc(3);       // Vector Fz
    Fx = gsl_vector_alloc(3);       // Vector Fx
    Fy = gsl_vector_alloc(3);       // Vector Fy
    acc = gsl_vector_alloc(3);      // Vector acc
    test_ug = gsl_matrix_alloc(9, 6001);
    test_xg = gsl_matrix_alloc(3, 6001);

}

/*
 * This Method returns the System vectors.
 * It gets an input int representing the vector:
 * 1  -> xg
 * 2  -> ug
 * 3  -> xg_alt
 * 4  -> ug_alt
 * 5  -> acc_alt
 * 6  -> delta_xg
 * 7  -> delta_ug
 * 8  -> alpha_r
 * 9  -> alpha_x
 * 10 -> alpha_y
 * 11 -> beta
 * 12 -> v
 * 13 -> vr
 * 14 -> sr
 * 15 -> mu
 * 16 -> mux
 * 17 -> muy
 * 18 -> acc
 */
gsl_vector* getVector(size_t n) {
    // Building the switch case structure

    switch (n) {
        case 1:
            return xg;
        case 2:
            return ug;
        case 3:
            return xg_alt;
        case 4:
            return ug_alt;
        case 5:
            return acc_alt;
        case 6:
            return delta_x;
        case 7:
            return delta_u;
        case 8:
            return alpha_r;
        case 9:
            return alpha_x;
        case 10:
            return alpha_y;
        case 11:
            return beta;
        case 12:
            return v;
        case 13:
            return vr;
        case 14:
            return sr;
        case 15:
            return mu;
        case 16:
            return mux;
        case 17:
            return muy;
        case 18:
            return acc;
        default:
            return 0;
    }
}

/*
 * This Method returns the System matrices.
 * It gets an input int representing the matrix:
 * 1 -> C
 * 2 -> D
 */
gsl_vector* getMatrix(size_t n){
    switch (n) {
        case 1:
            return C;
        case 2:
            return D;
        default:
            return 0;
    }
}

/*
 * For debugging purposes
 * All Vectors and matrices should be saved here
 */
void testVector(size_t n){

    gsl_matrix_get_col(ug, test_ug, n);
    gsl_matrix_get_col(xg, test_xg, n);
}

void initTest(){

    // Creating matrix test_ug for 60 sec simulation
    gsl_matrix_set_row(test_ug, 0, saving(1));      // col with 6001 values of n_1
    gsl_matrix_set_row(test_ug, 1, saving(2));      // col with 6001 values of n_2
    gsl_matrix_set_row(test_ug, 2, saving(3));      // col with 6001 values of n_3
    gsl_matrix_set_row(test_ug, 3, saving(7));      // col with 6001 values of delta_dyn_1
    gsl_matrix_set_row(test_ug, 4, saving(8));      // col with 6001 values of delta_dyn_2
    gsl_matrix_set_row(test_ug, 5, saving(9));      // col with 6001 values of delta_dyn_3
    gsl_matrix_set_row(test_ug, 6, saving(4));      // col with 6001 values of delta_1
    gsl_matrix_set_row(test_ug, 7, saving(5));      // col with 6001 values of delta_2
    gsl_matrix_set_row(test_ug, 8, saving(6));      // col with 6001 values of delta_3

    gsl_matrix_set_row(test_xg, 0, saving(10));     // col with 6001 values of v_x
    gsl_matrix_set_row(test_xg, 1, saving(11));     // col with 6001 values of v_y
    gsl_matrix_set_row(test_xg, 2, saving(12));     // col with 6001 values of psi_p
}

/*
 * This method calculates the vectors alpha_r, alpha_x and alpha_y and
 * saves it.
 * The formula is taken from Jan's Horizontal model line 25 to 27
 */
void slipage() {

    // Calculation of Vector indices alpha_r
    gsl_vector_set(alpha_r, 0, gsl_vector_get(ug, 3) - gsl_vector_get(ug, 6));      // Index 1
    gsl_vector_set(alpha_r, 1, gsl_vector_get(ug, 4) - gsl_vector_get(ug, 7));      // Index 2
    gsl_vector_set(alpha_r, 2, gsl_vector_get(ug, 5) - gsl_vector_get(ug, 8));      // Index 3

    // Calculation of Vector indices alpha_x
    gsl_vector_set(alpha_x, 0, gsl_vector_get(alpha_r, 0) * sin(gsl_vector_get(ug, 3)));       // Index 1
    gsl_vector_set(alpha_x, 1, gsl_vector_get(alpha_r, 1) * sin(gsl_vector_get(ug, 4)));       // Index 2
    gsl_vector_set(alpha_x, 2, gsl_vector_get(alpha_r, 2) * sin(gsl_vector_get(ug, 5)));       // Index 3

    // Calculation of Vector indices alpha_y
    gsl_vector_set(alpha_y, 0, gsl_vector_get(alpha_r, 0) * cos(gsl_vector_get(ug, 3)));       // Index 1
    gsl_vector_set(alpha_y, 1, gsl_vector_get(alpha_r, 1) * cos(gsl_vector_get(ug, 4)));       // Index 2
    gsl_vector_set(alpha_y, 2, gsl_vector_get(alpha_r, 2) * cos(gsl_vector_get(ug, 5)));       // Index 3

}

/*
 * This method calculates the velocity of the driving simulator given by ADMA sensor Data
 * Formula taken from Jan's Horizontal model in line 30 -40
 * Function get access on xg and do the calculations
 */
void adma_velocity() {

    // Current values of vector indices v_x, v_y and psi_p
    double v_x = gsl_vector_get(xg, 0);
    double v_y = gsl_vector_get(xg, 1);
    double psi_p = gsl_vector_get(xg, 2);

    // Calculation of velocity v1 and angle beta1
    v_1 = sqrt( pow(v_x ,2) + pow((v_y + psi_p* LAENGE / SQRT3),2 ));
    beta_1 = atan2((v_y + psi_p* LAENGE / SQRT3),v_x);

    // Calculation of velocity v2 and angle beta2
    v_2 = sqrt(pow((v_x + (psi_p * LAENGE / SQRT3) * cos(7 * M_PI / 6) ),2) + pow((v_y + (psi_p * LAENGE / SQRT3) * sin(7 * M_PI / 6) ),2));
    beta_2 = atan2((v_y + (psi_p * LAENGE / SQRT3) * sin(7 * M_PI / 6) ), (v_x + (psi_p * LAENGE / SQRT3) * cos(7 * M_PI / 6) ));

    // Calculation of velocity v3 and angle beta3
    v_3 = sqrt(pow((v_x + (psi_p * LAENGE / SQRT3) * cos(11 * M_PI / 6) ),2) + pow((v_y + (psi_p * LAENGE / SQRT3) * sin(11 * M_PI / 6) ),2));
    beta_3 = atan2((v_y + (psi_p * LAENGE / SQRT3) * sin(11 * M_PI / 6) ), (v_x + (psi_p * LAENGE / SQRT3) * cos(11 * M_PI / 6) ));

    // TODO Unbedingt ueberpruefen ob die Werte in Vector form vorliegen muessen
    // Calculation of beta vector and velocity vector
    // Setting the values in vector beta and v
    gsl_vector_set(beta, 0, beta_1);
    gsl_vector_set(beta, 1, beta_2);
    gsl_vector_set(beta, 2, beta_3);
    gsl_vector_set(v, 0, v_1* fabs(cos(beta_1)));
    gsl_vector_set(v, 1, v_2* fabs(cos(beta_2)));
    gsl_vector_set(v, 1, v_3* fabs(cos(beta_3)));
}

/*
 * This method gets the vector indices of ug and calculates the slip
 * Formula taken from Jan's horizontal model in line 46 - 54
 * Function get access on the values of vector ug and set into
 * vector vr after some Calculation
 * Two cases are calculated, one for Propulsion and one for braking
 */
void slip() {

    /* Velocity of tire calculated by its rotation speed
     * Get value of ug divide it by 30 and multiply by radius
     * of tire and Pi
     */
    gsl_vector_set(vr, 0, (M_PI * RADIUS * gsl_vector_get(ug, 0)/ 30));
    gsl_vector_set(vr, 1, (M_PI * RADIUS * gsl_vector_get(ug, 1)/ 30));
    gsl_vector_set(vr, 2, (M_PI * RADIUS * gsl_vector_get(ug, 2)/ 30));

    // TODO herausfinden was der Vector sr fuer eine Bedeutung hat
    for (size_t i = 0; i <3 ; i++) {
        if (gsl_vector_get(vr,i) >= gsl_vector_get(v,i))
            // Propulsion
            gsl_vector_set(sr, i, (1 - gsl_vector_get(v, i)/ gsl_vector_get(vr, i)) );
        else
            // Brake
            gsl_vector_set(sr, i, (1 - gsl_vector_get(vr, i)/ gsl_vector_get(v, i)) );
    }
}

/*
 * Method calculates the Friction for x and y-axes
 * Function gets access on values of vector sr and
 * saves information into vector mu and later on in
 * vector mux and muy
 */
void friction() {

    // Filling up vector mu with linearized slip KS and indices of vector sr
    gsl_vector_set(mu, 0, KMU * gsl_vector_get(sr, 0));
    gsl_vector_set(mu, 1, KMU * gsl_vector_get(sr, 1));
    gsl_vector_set(mu, 2, KMU * gsl_vector_get(sr, 2));

    // Filling the vector mx with indices of vector mu and cosine of vector ug
    gsl_vector_set(mux, 0, gsl_vector_get(mu,0) * cos(gsl_vector_get(ug, 3)));
    gsl_vector_set(mux, 1, gsl_vector_get(mu,1) * cos(gsl_vector_get(ug, 4)));
    gsl_vector_set(mux, 2, gsl_vector_get(mu,2) * cos(gsl_vector_get(ug, 5)));

    // Filling the vector my with indices of vector mu and sinus of vector ug
    gsl_vector_set(muy, 0, gsl_vector_get(mu,0) * sin(gsl_vector_get(ug, 3)));
    gsl_vector_set(muy, 1, gsl_vector_get(mu,1) * sin(gsl_vector_get(ug, 4)));
    gsl_vector_set(muy, 2, gsl_vector_get(mu,2) * sin(gsl_vector_get(ug, 5)));


}

/*
 * This method calculates the acceleration of in x-direction and returns double ax
 * It calculates the equitation of motion with the formula of Jan Steier
 * The goal of this Function is pursued by some auxiliary functions in it
 * One of these are the summarized doubles mux_all, muy_all, alpha_x_all and
 * alpha_y_all to get a better overview
 * Then the big equitation (size ~ length of A3 page) is divided into smaller
 * parts like sums, products and fraction with nominator and denominator for
 * better understanding
 * For more Information look up in Jan's MasterThesis or line 63 in Horizontal model
 */
// TODO alle Einträge nochmal uberpruefen
double Bewegungsgleichung_ax() {

    // Declaration of output acceleration ax
    double ax;

    /*
     * For an overview and summary of the long equitation some auxiliary scalars are calculated
     * and stored into the following doubles
     * mux_all, muy_all, alpha_x_all, alpha_y_all, firstPart,
     */
    double mux_all = gsl_vector_get(mux,0) + gsl_vector_get(mux, 1) + gsl_vector_get(mux, 2);
    double muy_all = gsl_vector_get(muy,0) + gsl_vector_get(muy, 1) + gsl_vector_get(muy, 2);
    double alphax_all = gsl_vector_get(alpha_x,0) + gsl_vector_get(alpha_x, 1) + gsl_vector_get(alpha_x, 2);
    double alphay_all = gsl_vector_get(alpha_y,0) + gsl_vector_get(alpha_y, 1) + gsl_vector_get(alpha_y, 2);

    /* The firstPart is a summary of the first to additions in the big equitation of motion
     * g/3 * [mu_x(0) + mu_x(1) + mu_x(2)] + ca / m * [alpha_x(0) + alpha_x(1) + alpha_x(2)]
     */
    double firstPart = GRAVY/3 * mux_all + C_a / MASSE * alphax_all;

    /* The Factor p1 is a summary of the first big factor in the big equitation
     * [hcg/l * (mux(2) - mux(3)] / [l / hcg  - muy(2) - muy(3)]
     */
    double p1 = (HCG/LAENGE * (gsl_vector_get(mux, 1) - gsl_vector_get(mux, 2)) ) / ( LAENGE / HCG - gsl_vector_get(muy, 2) + gsl_vector_get(muy, 2));

    /* The factor p2 is a summary of the second big factor in the fraction
     * mu_y(0) + mu_y(1) + mu_y(2) + ca / m * [alpha_y(0) + alpha_y(1) + alpha_y(2)
     */
    double p2 = muy_all + C_a / MASSE * alphay_all;

    // The numerator of the big fraction
    double fraction_1 = firstPart + p1 * p2;

    // The first quotient part for the big denominator
    double q1 = 1 - HCG / LAENGE / SQRT3 * (gsl_vector_get(mux, 1) + gsl_vector_get(mux, 2) - 2*gsl_vector_get(mux, 0)) - HCG/LAENGE * (gsl_vector_get(mux,1) - gsl_vector_get(mux,2));

    // The second quotient part for the big denominator
    double q2 = (LAENGE / HCG - gsl_vector_get(muy, 2) + gsl_vector_get(muy, 2)) * (HCG / LAENGE * (gsl_vector_get(muy, 1) + gsl_vector_get(muy, 2) - 2 * gsl_vector_get(muy, 0)) );

    // The denominator of the big fraction
    double fraction_2 = q1 / q2 ;

    ax = fraction_1 / fraction_2;

    return ax;
}


/*
 * This Method calculates the acceleration ay.
 * It describes the acceleration in y-axes and is important for the understanding
 * of vehicle motion. With ay you could reconstruct how the vehicle moves
 * and which action forces affect it.
 * It Solves the equitation with Jan's Formula from Horizontal model.
 * For getting the result, Bewegungsgleichung_ax() must be called
 */
double Bewegungsgleichung_ay() {

    // Declaration of double ay
    double ay;

    // Calculation of auxiliary doubles for further use in line below
    double muy_all = gsl_vector_get(muy,0) + gsl_vector_get(muy, 1) + gsl_vector_get(muy, 2);
    double alphay_all = gsl_vector_get(alpha_y,0) + gsl_vector_get(alpha_y, 1) + gsl_vector_get(alpha_y, 2);

    /*
     * ay is calculated by a big fraction
     * For better understanding look up the equitation in Jan's Master Thesis
     * or his Matlab Code in Horizontal Model
     */
    ay = (GRAVY/3 * muy_all + HCG / LAENGE / SQRT3 * (gsl_vector_get(muy,1) + gsl_vector_get(muy, 2) - 2 * gsl_vector_get(muy, 0))
            * Bewegungsgleichung_ax() + C_a / MASSE * alphay_all) / (1 - HCG / LAENGE *  (gsl_vector_get(muy, 1) - gsl_vector_get(muy,2)));
    return ay;
}

/*
 * This Method calculates the contact Forces of each point in the triangle
 * and returns a vector containing all three forces
 */
void AufstandsKraefte() {

    // current saving of ax and ay for caluclations
    double a_x = Bewegungsgleichung_ax();
    double a_y = Bewegungsgleichung_ay();
    /* Declaration of output vector and
     * Calculation of the three Forces FZ_1, FZ_2 and FZ_3
     */
    double FZ_1 = MASSE * GRAVY / 3 - 2 * MASSE * HCG / SQRT3 / LAENGE * a_x;
    double FZ_2 = MASSE * GRAVY / 3 + MASSE * HCG / LAENGE * (a_x / SQRT3 - a_y);
    double FZ_3 = MASSE * GRAVY / 3 + MASSE * HCG / LAENGE * (a_x / SQRT3 + a_y);

    // Putting the three forces into the vector
    gsl_vector_set(Fz,0, FZ_1);
    gsl_vector_set(Fz,1, FZ_2);
    gsl_vector_set(Fz,2, FZ_3);


}

/*
 * This Method calculates the Forces on the Tire and puts it into the Vector
 * Fx and Fy.
 * The formula is taken from Jan's Horizontal model
 * Fx(x) = Fz(x) * mux(x) + C_a * alpha_x(x)
 */
void RadKraefte() {
    // Putting the calculated scalar into the vector Fx
    gsl_vector_set(Fx,0, gsl_vector_get(Fz,0) * gsl_vector_get(mux,0) + C_a * gsl_vector_get(alpha_x,0));
    gsl_vector_set(Fx,1, gsl_vector_get(Fz,1) * gsl_vector_get(mux,1) + C_a * gsl_vector_get(alpha_x,1));
    gsl_vector_set(Fx,2, gsl_vector_get(Fz,2) * gsl_vector_get(mux,2) + C_a * gsl_vector_get(alpha_x,2));

    // Putting the calculated scalar into the vector Fy
    gsl_vector_set(Fy,0, gsl_vector_get(Fz,0) * gsl_vector_get(muy,0) + C_a * gsl_vector_get(alpha_y,0));
    gsl_vector_set(Fy,1, gsl_vector_get(Fz,1) * gsl_vector_get(muy,1) + C_a * gsl_vector_get(alpha_y,1));
    gsl_vector_set(Fy,2, gsl_vector_get(Fz,2) * gsl_vector_get(muy,2) + C_a * gsl_vector_get(alpha_y,2));
}

/*
 * This Method calculate the yaw acceleration and returns the general acceleration
 * in Form of the vector acc
 */
void GierbewegungBerechnen() {

    // Calculation of psi_pp with formula taken from Jan's Modell
    double psi_pp = LAENGE/THETA * (gsl_vector_get(Fy,0) / SQRT3 - gsl_vector_get(Fy,1) / 2 / SQRT3 - gsl_vector_get(Fy,2)
            / 2 / SQRT3 - gsl_vector_get(Fx, 1) / 2 + gsl_vector_get(Fx,2) / 2);

    // Filling the Vector acc with acceleration ax, ay and psi_pp
    gsl_vector_set(acc, 0, Bewegungsgleichung_ax());
    gsl_vector_set(acc, 1, Bewegungsgleichung_ay());
    gsl_vector_set(acc, 2, psi_pp);
    /*
    printf("  acc--> ");
    printf("%f ", gsl_vector_get(acc,0));
    printf("%f ", gsl_vector_get(acc,1));
    printf("%f ", gsl_vector_get(acc,2));
    printf("  <--- acc \n");
    */
}

/*
 * This Method calculates the Big System matrices C and D.
 * For the equitation some auxiliary parameters are calculated.
 *
 */
void SystemmatrixBerechnen() {

    // Calculating auxiliary parameters acc for further calculations
    double acc_one = gsl_vector_get(acc, 0) - gsl_vector_get(acc_alt, 0);
    double acc_two = gsl_vector_get(acc, 1) - gsl_vector_get(acc_alt, 1);
    double acc_three = gsl_vector_get(acc, 2) - gsl_vector_get(acc_alt, 2);

    // Calculating auxiliary parameters xg for further calculations
    double xg_one = gsl_vector_get(xg, 0) - gsl_vector_get(xg_alt, 0);
    double xg_two = gsl_vector_get(xg, 1) - gsl_vector_get(xg_alt, 1);
    double xg_three = gsl_vector_get(xg, 2) - gsl_vector_get(xg_alt, 2);

    // Calculating auxiliary parameters ug for further calculations
    double ug_one = gsl_vector_get(ug, 0) - gsl_vector_get(ug_alt, 0);
    double ug_two = gsl_vector_get(ug, 1) - gsl_vector_get(ug_alt, 1);
    double ug_three = gsl_vector_get(ug, 2) - gsl_vector_get(ug_alt, 2);
    double ug_four = gsl_vector_get(ug, 3) - gsl_vector_get(ug_alt, 3);
    double ug_five = gsl_vector_get(ug, 4) - gsl_vector_get(ug_alt, 4);
    double ug_six = gsl_vector_get(ug, 5) - gsl_vector_get(ug_alt, 5);
    double ug_seven = gsl_vector_get(ug, 6) - gsl_vector_get(ug_alt, 6);
    double ug_eight = gsl_vector_get(ug, 7) - gsl_vector_get(ug_alt, 7);
    double ug_nine = gsl_vector_get(ug, 8) - gsl_vector_get(ug_alt, 8);

    //printf("%f, %f , %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n", acc_one, acc_two, acc_three, xg_one, xg_two, xg_three, ug_one, ug_two, ug_three, ug_four, ug_five, ug_six, ug_seven, ug_eight, ug_nine);

    // TODO Unbedingt kommentieren
    /*
     * This for loop calculates the three rows in the Matrix D
     * In the first if-Statement the difference between xg and xg_alt
     * leads to an increase of xg_alt if the difference is not zero
     * After updating xg_alt the values of the D matrix are calculated
     * and put into special places.
     * The formula for each indices is written in matrix_setter
     * For further information look in manual gsl_matrix_set or
     * Jan's Horizontal model
     */
    for (size_t i = 0; i <3 ; i++) {

        // Proof if difference between xg and xg_alt is zero or not and increase xg_alt by 0.0001
        if ((gsl_vector_get(xg, i) - gsl_vector_get(xg_alt, i)) != 0)
            gsl_vector_set(xg_alt, i, gsl_vector_get(xg_alt, i) + 0.0001);

        // Switch case puts fraction of acc and xg into the matrix
        switch (i) {
            // for i equals 0 fill 4th, 10th and 16th index of D
            case 0:
                gsl_vector_set(D, 3, acc_one / xg_one);
                gsl_vector_set(D, 9, acc_two / xg_one);
                gsl_vector_set(D, 15, acc_three / xg_one);
            // for i equals 1 fill 5th, 11th and 17th index of D
            case 1:
                gsl_vector_set(D, 4, acc_one / xg_two);
                gsl_vector_set(D, 10, acc_two / xg_two);
                gsl_vector_set(D, 16, acc_three / xg_two);
            // for i equals 1 fill 6th, 12th and 18th index of D
            case 2:
                gsl_vector_set(D, 5, acc_one / xg_three);
                gsl_vector_set(D, 11, acc_two / xg_three);
                gsl_vector_set(D, 17, acc_three / xg_three);
            default:
                break;
        }
    }

    /*
     * This for loop does the same as the loop above with some differences.
     * The size of the counter is 3 times bigger then before and here
     * it also calculates indices of botn Matrices, C and D
     * Therefore the switch-case block is bigger then before and
     * the auxiliary vectors ug are more then xg
     */
    for (size_t j = 0; j <9 ; j++) {

        // Proof if difference between ug and ug_alt is 0
        if ((gsl_vector_get(ug,j) - gsl_vector_get(ug_alt,j)) != 0)
            // if true increase value of index j of ug_alt by 0.0001
            gsl_vector_set(ug_alt, j, gsl_vector_get(ug_alt,j) + 0.0001);
        // switch separate the calculation depending on counter j
        switch (j){

            // if j is 0 set those values in matrix C
            case 0:
                gsl_vector_set(C, 3, acc_one / ug_one);
                gsl_vector_set(C, 9, acc_two / ug_one);
                gsl_vector_set(C, 15, acc_three / ug_one);
            // if j is 1 set those values in matrix C
            case 1:
                gsl_vector_set(C, 4, acc_one / ug_two);
                gsl_vector_set(C, 10, acc_two / ug_two);
                gsl_vector_set(C, 16, acc_three / ug_two);
            // if j is 2 set those values in matrix C
            case 2:
                gsl_vector_set(C, 5, acc_one / ug_three);
                gsl_vector_set(C, 11, acc_two / ug_three);
                gsl_vector_set(C, 17, acc_three / ug_three);
            // if j is 3 set those values in matrix C
            case 3:
                gsl_vector_set(C, 0, acc_one / ug_four);
                gsl_vector_set(C, 6, acc_two / ug_four);
                gsl_vector_set(C, 12, acc_three / ug_four);
            // if j is 4 set those values in matrix C
            case 4:
                gsl_vector_set(C, 1, acc_one / ug_five);
                gsl_vector_set(C, 7, acc_two / ug_five);
                gsl_vector_set(C, 13, acc_three / ug_five);
            // if j is 5 set those values in matrix C
            case 5:
                gsl_vector_set(C, 2, acc_one / ug_six);
                gsl_vector_set(C, 8, acc_two / ug_six);
                gsl_vector_set(C, 14, acc_three / ug_six);
            // if j is 6 set those values in matrix D
            case 6:
                gsl_vector_set(D, 0, acc_one / ug_seven);
                gsl_vector_set(D, 6, acc_two / ug_seven);
                gsl_vector_set(D, 12, acc_three / ug_seven);
            // if j is 7 set those values in matrix D
            case 7:
                gsl_vector_set(D, 1, acc_one / ug_eight);
                gsl_vector_set(D, 7, acc_two / ug_eight);
                gsl_vector_set(D, 13, acc_three / ug_eight);
            case 8:
                gsl_vector_set(D, 2, acc_one / ug_nine);
                gsl_vector_set(D, 8, acc_two / ug_nine);
                gsl_vector_set(D, 14, acc_three / ug_nine);
            default:
                break;
        }
    }

    /*
    printf("acc_alt_1  |  acc_alt_2  |  acc_alt_3\n");
    printf(" %f |", gsl_vector_get(acc_alt,0));
    printf(" %f   |", gsl_vector_get(acc_alt,1));
    printf(" %f \n", gsl_vector_get(acc_alt,2));
    */
    // Set acc as the new acc_alt
    gsl_vector_memcpy(acc_alt, acc);

    /*
    printf("acc_1      |     acc_2   |   acc_3\n");
    printf(" %f |", gsl_vector_get(acc,0));
    printf(" %f   |", gsl_vector_get(acc,1));
    printf(" %f \n", gsl_vector_get(acc,2));
     */
}


/*
 * Calculation of the Delta Vectors
 * Storing the current value of xg and ug in delta_xg
 * and delta_ug.
 * Method subtract xg_alt/ug_alt from deltas and saves
 * it in the delta vectors
 */
void deltasBerechnen(){
    delta_x = xg;       // Copy of vector xg
    delta_u = ug;       // Copy of vector ug
    gsl_vector_sub(delta_x, xg_alt);
    gsl_vector_sub(delta_u, ug_alt);
}


/*
 * This method copies the current values of ug and xg into
 * ug_alt and xg_alt
 */
void saving_current_state(){
    gsl_vector_memcpy(xg_alt, xg);
    gsl_vector_memcpy(ug_alt,ug);
    gsl_vector_memcpy(acc_alt,acc);
}