//
// Created by Farid Bonakdar on 30.06.16.
//

#include "horizontalModel.h"
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


#define KS 5        // Linearisierung mu-Schlupf war 5
#define C_a 9095    // Schräglaufsteifigkeit
#define R 0.15      // Radradius
#define M 1200      // total mass kg ursprünglich 308.52
#define HCG 0.48303 // height of the center of gravity m  -- CarMaker// 0.55m platform //total 0.48303 selbt gerechnet mit allen Komponenten aus CarMaker
#define L 2.28      // equilateral triangle setup m ursprünglich 2.54
#define G 9.81      // gravity constant
#define THETA 944.8465962
#define SQRT3 1.73205

/*
 * Declaration of Systemmatrices and inner Systemvector for further use in methods
 */

    gsl_vector* xg;       // Vector xg
    gsl_vector* ug;       // Vector ug
    gsl_vector* xg_alt;   // Vector xg_alt
    gsl_vector* ug_alt;   // Vector ug_Alt
    gsl_vector* acc_alt;  // Vector acc_alt
    gsl_vector* delta_x;  // Vector delta_x
    gsl_vector* delta_u;  // Vector delta_u
    gsl_vector* alpha_r;  // Vector alpha_r
    gsl_vector* alphax;   // Vector alphax
    gsl_vector* alphay;   // Vector alphay
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
    gsl_vector* mux;       // Vector mu
    gsl_vector* muy;       // Vector mu


/*
 * Initialization of Systemmatrix C and D.
 * Allocation of inner Systemvector xg, ug, xg_alt, ug_alt and acc_alt
 * The size of vector and matrix taken from Jans Matlab Code (Horizontalmodell)
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
    alphax = gsl_vector_alloc(3);   // Vector alphax
    alphay = gsl_vector_alloc(3);   // Vector alphay
    C = gsl_matrix_calloc(18,1);    // Matrix C
    D = gsl_matrix_calloc(18,1);    // Matrix D
    beta = gsl_vector_alloc(3);     // Vector beta
    v = gsl_vector_alloc(3);        // Vector v
    vr = gsl_vector_alloc(3);       // Vector vr
    sr = gsl_vector_alloc(3);       // Vector sr
    mu = gsl_vector_alloc(3);       // Vector mu
    mux = gsl_vector_alloc(3);       // Vector mx
    muy = gsl_vector_alloc(3);       // Vector my

}

/*
 * This Method returns the Systemvectors.
 * It gets an input int representing the vector:
 * 1  -> xg
 * 2  -> ug
 * 3  -> xg_alt
 * 4  -> ug_alt
 * 5  -> acc_alt
 * 6  -> delta_xg
 * 7  -> delta_ug
 * 8  -> alpha_r
 * 9  -> alphax
 * 10 -> alphay
 */
gsl_vector* getVector(int n) {
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
            return alphax;
        case 10:
            return alphay;
        default:
            break;
    }
}

/*
 * This Method returns the Systemmatrices.
 * It gets an input int representing the matrix:
 * 1 -> C
 * 2 -> D
 */
gsl_matrix* getMatrix(int n){
    switch (n) {
        case 1:
            return C;
        case 2:
            return D;
        default:
            break;
    }
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
 * This method calculates the vectors alpha_r, alphax and alphay and
 * saves it.
 * The formula is taken from Jans Horizontalmodell line 25 to 27
 */
void schraeglaufwinkel() {

    // Calculation of Vector indices alpha_r
    gsl_vector_set(alpha_r, 0, gsl_vector_get(ug, 3) - gsl_vector_get(ug, 6));      // Index 1
    gsl_vector_set(alpha_r, 1, gsl_vector_get(ug, 4) - gsl_vector_get(ug, 7));      // Index 2
    gsl_vector_set(alpha_r, 2, gsl_vector_get(ug, 5) - gsl_vector_get(ug, 8));      // Index 3

    // Calculation of Vector indices alphax
    gsl_vector_set(alphax, 0, gsl_vector_get(alpha_r, 0) * sin(gsl_vector_get(ug, 3)));       // Index 1
    gsl_vector_set(alphax, 1, gsl_vector_get(alpha_r, 1) * sin(gsl_vector_get(ug, 4)));       // Index 2
    gsl_vector_set(alphax, 2, gsl_vector_get(alpha_r, 2) * sin(gsl_vector_get(ug, 5)));       // Index 3

    // Calculation of Vector indices alphay
    gsl_vector_set(alphay, 0, gsl_vector_get(alpha_r, 0) * cos(gsl_vector_get(ug, 3)));       // Index 1
    gsl_vector_set(alphay, 1, gsl_vector_get(alpha_r, 1) * cos(gsl_vector_get(ug, 4)));       // Index 2
    gsl_vector_set(alphay, 2, gsl_vector_get(alpha_r, 2) * cos(gsl_vector_get(ug, 5)));       // Index 3

}

/*
 * This method calculates the velocity of the driving simulator given by ADMA sensor Data
 * Formula taken from Jans Horizontalmodell in line 30 -40
 */
void Adma_geschwindigkeit() {

    // Current vector indices v_x, v_y and psi_p
    double v_x = gsl_vector_get(xg, 0);
    double v_y = gsl_vector_get(xg, 1);
    double psi_p = gsl_vector_get(xg, 2);

    // Calculation of velocity v1 and angle beta1
    v_1 = sqrt( pow(v_x ,2) + pow((v_y + psi_p* L / SQRT3),2 ));
    beta_1 = atan2((v_y + psi_p* L / SQRT3),v_x);

    // Calculation of velcocity v2 and angle beta2
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
 * Formula taken from Jans horizontalmodell in line 46 - 54
 */
void schlupfBerechnung() {

    // Velocity of tire calculatet by its rotation speed
    gsl_vector_set(vr, 0, (M_PI * R * gsl_vector_get(ug, 0)/ 30));
    gsl_vector_set(vr, 1, (M_PI * R * gsl_vector_get(ug, 1)/ 30));
    gsl_vector_set(vr, 2, (M_PI * R * gsl_vector_get(ug, 2)/ 30));

    // TODO herausfinden was der Vector sr fuer eine Bedeutung hat
    for (int i = 0; i <3 ; i++) {
        if (gsl_vector_get(vr,i) >= gsl_vector_get(v,i))
            // Antriebsfall
            gsl_vector_set(sr, i, (1 - gsl_vector_get(v, i)/ gsl_vector_get(vr, i)) );
        else
            // Bremsfall
            gsl_vector_set(sr, i, (1 - gsl_vector_get(vr, i)/ gsl_vector_get(v, i)) );
    }
}

/*
 * Method calculates the Friction for x and y-axes
 */
void ReibwertBerechnung() {

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
 * For more Information look up in Jan's MasterThesis or line 63 in Horizontalmodell
 */
// TODO alle Einträge nochmal uberpruefen
double Bewegungsgleichung_ax() {
    // Declaration of output acceleration ax
    double ax;
    /*
     * For an overview and summary of the long equitation some auxiliary scalars are calculated
     * and stored into the following doubles
     * mux_all, muy_all, alphax_all, alphay_all, firstPart,
     * p1 and p2
     *
     */
    double mux_all = gsl_vector_get(mux,0) + gsl_vector_get(mux, 1) + gsl_vector_get(mux, 2);
    double muy_all = gsl_vector_get(muy,0) + gsl_vector_get(muy, 1) + gsl_vector_get(muy, 2);
    double alphax_all = gsl_vector_get(alphax,0) + gsl_vector_get(alphax, 1) + gsl_vector_get(alphax, 2);
    double alphay_all = gsl_vector_get(alphay,0) + gsl_vector_get(alphay, 1) + gsl_vector_get(alphay, 2);

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


double Bewegungsgleichung_ay() {
    double ay;

    return ay;
}

double AufstandsKraefte() {
    return 0;
}

double RadKraefte() {
    return 0;
}

double GierbewegungBerechnen() {
    return 0;
}

double SystemmatrixBerechnen() {
    return 0;
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
