//
// Created by Farid Bonakdar on 30.06.16.
//
#include "linear.h"

#define SIGNALE 42.0
#define REGLER 'zahl'
#define MAT "arr"
#define USW "NOCHMEHR"


// TODO Funktion rausschmeissen oder ersetzen
double fun(double one, double two){
    return SIGNALE;
}

//TODO Methode zur Berechnung der K_p Matrix
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

