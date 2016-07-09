//
// Created by Farid Bonakdar on 30.06.16.
//

#include "horizontalModel.h"
#include <math.h>


#define KS 5                        // linearized mu-slip about 5
#define C_a 9095                    // slippage
#define R 0.15                      // radius of tire
#define M 1200                      // total mass kg originally 308.52
#define HCG 0.48303                 // height of the center of gravity m  -- CarMaker self calculated // 0.55m platform
#define L 2.28                      // equilateral triangle setup m initially 2.54
#define G 9.81                      // gravity constant
#define THETA 944.8465962           // Angle Theta
#define SQRT3 1.73205               // square root of 3

/*
 * Declaration of System matrices and inner System vector for further use in methods.
 * Stored here to get acces to each scalar, vector and matrix.
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
    gsl_matrix* C;        // Matrix C
    gsl_matrix* D;        // Matrix D
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
    C = gsl_matrix_calloc(18,1);    // Matrix C
    D = gsl_matrix_calloc(18,1);    // Matrix D
    beta = gsl_vector_alloc(3);     // Vector beta
    v = gsl_vector_alloc(3);        // Vector v
    vr = gsl_vector_alloc(3);       // Vector vr
    sr = gsl_vector_alloc(3);       // Vector sr
    mu = gsl_vector_alloc(3);       // Vector mu
    mux = gsl_vector_alloc(3);      // Vector mx
    muy = gsl_vector_alloc(3);      // Vector my
    Fz = gsl_vector_alloc(3);       // Vector Fz
    Fx = gsl_vector_alloc(3);       // Vector Fx
    Fy = gsl_vector_alloc(3);       // Vector Fy
    acc = gsl_vector_alloc(3);      // Vector acc

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
 */
gsl_vector* getVector(int n) {
    // Declaration of output vector
    gsl_vector* out;

    // If vector ug, ug_alt or delta_u is called allocate
    // out with size 9
    // Else allocate out with size 3
    if (n == 2 || n == 4 || n == 7)
        out = gsl_vector_alloc(9);
    else
        out = gsl_vector_alloc(3);

    switch (n) {
        case 1:
            out = xg;
        case 2:
            out = ug;
        case 3:
            out = xg_alt;
        case 4:
            out = ug_alt;
        case 5:
            out = acc_alt;
        case 6:
            out = delta_x;
        case 7:
            out = delta_u;
        case 8:
            out = alpha_r;
        case 9:
            out = alpha_x;
        case 10:
            out = alpha_y;
        case 11:
            out = beta;
        case 12:
            out = v;
        case 13:
            out = vr;
        case 14:
            out = sr;
        case 15:
            out = mu;
        case 16:
            out = mux;
        case 17:
            out = muy;
        default:
            break;
    }
    return out;
}

/*
 * This Method returns the System matrices.
 * It gets an input int representing the matrix:
 * 1 -> C
 * 2 -> D
 */
gsl_matrix* getMatrix(int n){
    gsl_matrix* out = gsl_matrix_alloc(18,1);
    switch (n) {
        case 1:
            out = C;
        case 2:
            out = D;
        default:
            break;
    }
    return out;
}

/*
 * For debugging purposes
 * All Vectors and matrices should be saved here
 */
void testVector(){

    for (double j = 0; j < 3 ; j++) {
        gsl_vector_set(xg, j, j + 2.3);
        gsl_vector_set(xg_alt, j, j + 2.3);
    }

    for (double i = 0; i < 9 ; i++) {
        gsl_vector_set(ug, i, i+2.7);
        gsl_vector_set(ug_alt, i, i + 3.5);
    }

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
 */
void adma_velocity() {

    // Current vector indices v_x, v_y and psi_p
    double v_x = gsl_vector_get(xg, 0);
    double v_y = gsl_vector_get(xg, 1);
    double psi_p = gsl_vector_get(xg, 2);

    // Calculation of velocity v1 and angle beta1
    v_1 = sqrt( pow(v_x ,2) + pow((v_y + psi_p* L / SQRT3),2 ));
    beta_1 = atan2((v_y + psi_p* L / SQRT3),v_x);

    // Calculation of velocity v2 and angle beta2
    v_2 = sqrt(pow((v_x + (psi_p * L / SQRT3) * cos(7 * M_PI / 6) ),2) + pow((v_y + (psi_p * L / SQRT3) * sin(7 * M_PI / 6) ),2));
    beta_2 = atan2((v_y + (psi_p * L / SQRT3) * sin(7 * M_PI / 6) ), (v_x + (psi_p * L / SQRT3) * cos(7 * M_PI / 6) ));

    // Calculation of velocity v3 and angle beta3
    v_3 = sqrt(pow((v_x + (psi_p * L / SQRT3) * cos(11 * M_PI / 6) ),2) + pow((v_y + (psi_p * L / SQRT3) * sin(11 * M_PI / 6) ),2));
    beta_3 = atan2((v_y + (psi_p * L / SQRT3) * sin(11 * M_PI / 6) ), (v_x + (psi_p * L / SQRT3) * cos(11 * M_PI / 6) ));

    // TODO Unbedingt ueberpruefen ob die Werte in Vector form vorliegen muessen
    // Calculation of beta vector and velocity vector
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
 */
void slip() {

    // Velocity of tire calculated by its rotation speed
    gsl_vector_set(vr, 0, (M_PI * R * gsl_vector_get(ug, 0)/ 30));
    gsl_vector_set(vr, 1, (M_PI * R * gsl_vector_get(ug, 1)/ 30));
    gsl_vector_set(vr, 2, (M_PI * R * gsl_vector_get(ug, 2)/ 30));

    // TODO herausfinden was der Vector sr fuer eine Bedeutung hat
    for (int i = 0; i <3 ; i++) {
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
 */
void friction() {

    // Filling up vector mu with linearized slip KS and indices of vector sr
    gsl_vector_set(mu, 0, KS * gsl_vector_get(sr, 0));
    gsl_vector_set(mu, 1, KS * gsl_vector_get(sr, 1));
    gsl_vector_set(mu, 2, KS * gsl_vector_get(sr, 2));

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
     * p1 and p2
     *
     */
    double mux_all = gsl_vector_get(mux,0) + gsl_vector_get(mux, 1) + gsl_vector_get(mux, 2);
    double muy_all = gsl_vector_get(muy,0) + gsl_vector_get(muy, 1) + gsl_vector_get(muy, 2);
    double alphax_all = gsl_vector_get(alpha_x,0) + gsl_vector_get(alpha_x, 1) + gsl_vector_get(alpha_x, 2);
    double alphay_all = gsl_vector_get(alpha_y,0) + gsl_vector_get(alpha_y, 1) + gsl_vector_get(alpha_y, 2);

    // The firstPart is a summary of the first to additions in the big equitation of motion
    // g/3 * [mu_x(0) + mu_x(1) + mu_x(2)] + ca / m * [alpha_x(0) + alpha_x(1) + alpha_x(2)]
    double firstPart = G/3 * mux_all + C_a / M * alphax_all;

    // The Factor p1 is a summary of the first big factor in the big equitation
    // [hcg/l * (mux(2) - mux(3)] / [l / hcg  - muy(2) - muy(3)]
    double p1 = (HCG/L * (gsl_vector_get(mux, 1) - gsl_vector_get(mux, 2)) ) / ( L / HCG - gsl_vector_get(muy, 2) + gsl_vector_get(muy, 2));

    // The factor p2 is a summary of the second big factor in the fraction
    // mu_y(0) + mu_y(1) + mu_y(2) + ca / m * [alpha_y(0) + alpha_y(1) + alpha_y(2)
    double p2 = muy_all + C_a / M * alphay_all;

    // The numerator of the big fraction
    double frac1 = firstPart + p1 * p2;

    // The first quotient part for the big denominator
    double q1 = 1 - HCG / L / SQRT3 * (gsl_vector_get(mux, 1) + gsl_vector_get(mux, 2) - 2*gsl_vector_get(mux, 0)) - HCG/L * (gsl_vector_get(mux,1) - gsl_vector_get(mux,2));

    // The second quotient part for the big denominator
    double q2 = (L / HCG - gsl_vector_get(muy, 2) + gsl_vector_get(muy, 2)) * (HCG / L * (gsl_vector_get(muy, 1) + gsl_vector_get(muy, 2) - 2 * gsl_vector_get(muy, 0)) );

    // The denominator of the big fraction
    double frac2 = q1 / q2 ;

    ax = frac1 / frac2;

    return ax;
}


/*
 * This Method calculates the acceleration ay.
 * It Solves the equitation with Jan's Formula from Horizontal model.
 * For getting the result, Bewegungsgleichung_ax() must be called
 */
double Bewegungsgleichung_ay() {
    double ay;
    double muy_all = gsl_vector_get(muy,0) + gsl_vector_get(muy, 1) + gsl_vector_get(muy, 2);
    double alphay_all = gsl_vector_get(alpha_y,0) + gsl_vector_get(alpha_y, 1) + gsl_vector_get(alpha_y, 2);

    ay = (G/3 * muy_all + HCG / L / SQRT3 * (gsl_vector_get(muy,1) + gsl_vector_get(muy, 2) - 2 * gsl_vector_get(muy, 0))
            * Bewegungsgleichung_ax() + C_a / M * alphay_all) / (1 - HCG / L *  (gsl_vector_get(muy, 1) - gsl_vector_get(muy,2)));
    return ay;
}

/*
 * This Method calculates the Forces of each point in the triangle
 * and returns a vector containing all three forces
 */
gsl_vector * AufstandsKraefte() {
    // Declaration of output vector
    // Calculation of the three Forces FZ_1, FZ_2 and FZ_3
    double FZ_1 = M * G / 3 - 2 * M * HCG / SQRT3 / L * Bewegungsgleichung_ax();
    double FZ_2 = M * G / 3 + M * HCG / L * (Bewegungsgleichung_ax() / SQRT3 - Bewegungsgleichung_ay());
    double FZ_3 = M * G / 3 + M * HCG / L * (Bewegungsgleichung_ax() / SQRT3 + Bewegungsgleichung_ay());

    // Putting the three forces into the vector
    gsl_vector_set(Fz,0, FZ_1);
    gsl_vector_set(Fz,1, FZ_2);
    gsl_vector_set(Fz,2, FZ_3);

    return Fz;
}

/*
 * This Method calculates the Forces on the Tire and puts it into the Vector
 * Fx and Fy.
 * The formula is taken from Jan's Horizontal model
 * Fx(x) = Fz(x) * mux(x) + C_a * alpha_x(x)
 */
void RadKraefte() {
    // Putting the calculated scalar into the vector Fx
    gsl_vector_set(Fx,0, gsl_vector_get(AufstandsKraefte(),0) * gsl_vector_get(mux,0) + C_a * gsl_vector_get(alpha_x,0));
    gsl_vector_set(Fx,1, gsl_vector_get(AufstandsKraefte(),1) * gsl_vector_get(mux,1) + C_a * gsl_vector_get(alpha_x,1));
    gsl_vector_set(Fx,2, gsl_vector_get(AufstandsKraefte(),2) * gsl_vector_get(mux,2) + C_a * gsl_vector_get(alpha_x,2));

    // Putting the calculated scalar into the vector Fy
    gsl_vector_set(Fy,0, gsl_vector_get(AufstandsKraefte(),0) * gsl_vector_get(muy,0) + C_a * gsl_vector_get(alpha_y,0));
    gsl_vector_set(Fy,1, gsl_vector_get(AufstandsKraefte(),1) * gsl_vector_get(muy,1) + C_a * gsl_vector_get(alpha_y,1));
    gsl_vector_set(Fy,2, gsl_vector_get(AufstandsKraefte(),2) * gsl_vector_get(muy,2) + C_a * gsl_vector_get(alpha_y,2));
}

/*
 * This Method calculate the yaw acceleration and returns the general acceleration
 * in Form of the vector acc
 */
void GierbewegungBerechnen() {

    // Calculation of psi_pp with formula taken from Jan's Modell
    double psi_pp = L/THETA * (gsl_vector_get(Fy,0) / SQRT3 - gsl_vector_get(Fy,1) / 2 / SQRT3 - gsl_vector_get(Fy,2)
            / 2 / SQRT3 - gsl_vector_get(Fx, 1) / 2 + gsl_vector_get(Fx,2) / 2);

    // Filling the Vector acc with acceleration ax, ay and psi_pp
    gsl_vector_set(acc, 0, Bewegungsgleichung_ax());
    gsl_vector_set(acc, 1, Bewegungsgleichung_ay());
    gsl_vector_set(acc, 2, psi_pp);

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
    for (int i = 0; i <3 ; i++) {

        // Proof if difference between xg and xg_alt is zero or not and increase xg_alt by 0.0001
        if ((gsl_vector_get(xg, i) - gsl_vector_get(xg_alt, i)) != 0)
            gsl_vector_set(xg_alt, i, gsl_vector_get(xg_alt, i) + 0.0001);

        // Switch case puts fraction of acc and xg into the matrix
        switch (i) {
            // for i equals 0 fill 4th, 10th and 16th index of D
            case 0:
                gsl_matrix_set(D, 3, 0, acc_one / xg_one);
                gsl_matrix_set(D, 9, 0, acc_two / xg_one);
                gsl_matrix_set(D, 15, 0, acc_three / xg_one);
            // for i equals 1 fill 5th, 11th and 17th index of D
            case 1:
                gsl_matrix_set(D, 4, 0, acc_one / xg_two);
                gsl_matrix_set(D, 10, 0, acc_two / xg_two);
                gsl_matrix_set(D, 16, 0, acc_three / xg_two);
            // for i equals 1 fill 6th, 12th and 18th index of D
            case 2:
                gsl_matrix_set(D, 5, 0, acc_one / xg_three);
                gsl_matrix_set(D, 11, 0, acc_two / xg_three);
                gsl_matrix_set(D, 17, 0, acc_three / xg_three);
            default:
                break;
        }
    }

    for (int j = 0; j <8 ; j++) {
        if ((gsl_vector_get(ug,j) - gsl_vector_get(ug_alt,j)) != 0)
            gsl_vector_set(ug_alt, j, gsl_vector_get(ug_alt,j) + 0.0001);
        switch (j){
            case 0:
                gsl_matrix_set(C,3,0, acc_one / ug_one);
                gsl_matrix_set(C,9,0, acc_two / ug_one);
                gsl_matrix_set(C,15,0, acc_three / ug_one);
            case 1:
                gsl_matrix_set(C,4,0, acc_one / ug_two);
                gsl_matrix_set(C,10,0, acc_two / ug_two);
                gsl_matrix_set(C,16,0, acc_three / ug_two);
            case 2:
                gsl_matrix_set(C,5,0, acc_one / ug_three);
                gsl_matrix_set(C,11,0, acc_two / ug_three);
                gsl_matrix_set(C,17,0, acc_three / ug_three);
            case 3:
                gsl_matrix_set(C,0,0, acc_one / ug_four);
                gsl_matrix_set(C,6,0, acc_two / ug_four);
                gsl_matrix_set(C,12,0, acc_three / ug_four);
            case 4:
                gsl_matrix_set(C,1,0, acc_one / ug_five);
                gsl_matrix_set(C,7,0, acc_two / ug_five);
                gsl_matrix_set(C,13,0, acc_three / ug_five);
            case 5:
                gsl_matrix_set(C,2,0, acc_one / ug_six);
                gsl_matrix_set(C,8,0, acc_two / ug_six);
                gsl_matrix_set(C,14,0, acc_three / ug_six);
            case 6:
                gsl_matrix_set(D,0,0, acc_one / ug_seven);
                gsl_matrix_set(D,6,0, acc_two / ug_seven);
                gsl_matrix_set(D,12,0, acc_three / ug_seven);
            case 7:
                gsl_matrix_set(C,1,0, acc_one / ug_eight);
                gsl_matrix_set(C,7,0, acc_two / ug_eight);
                gsl_matrix_set(C,13,0, acc_three / ug_eight);
            case 8:
                gsl_matrix_set(C,2,0, acc_one / ug_nine);
                gsl_matrix_set(C,8,0, acc_two / ug_nine);
                gsl_matrix_set(C,14,0, acc_three / ug_nine);
            default:
                break;
        }
    }
    acc_alt = acc;
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
