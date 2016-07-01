#include <stdio.h>
#include "modelAntrieb.h"
#include "DGL_Berechnen.h"
#include "Subsysteme.h"
#include "calcMatrix_A_B.h"
#include "horizontalModel.h"
#include "linear.h"
#include "Lenkmotoren.h"
#include <time.h>



int main() {
    clock_t begin = clock();
    //double a = 5.0;
    printf("Auch hallo\n");
    printf("%f\n" , fun(2.0,3.0));
    // TODO funktioniert nicht !!!
    //printf("%4s", ausgabe());
    clock_t  end = clock();
    double time = (double) 1000*(end-begin)/CLOCKS_PER_SEC;
    printf("Berechnung dauert %.3f ms\n", time);
    return 0;
}