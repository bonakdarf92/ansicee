cmake_minimum_required(VERSION 3.5)
project(ansicee)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/include/gsl -lgsl")

include_directories(-I/usr/include/gsl)
set(SOURCE_FILES main.c linear.h linear.c modelAntrieb.c modelAntrieb.h calcMatrix_A_B.c calcMatrix_A_B.h horizontalModel.c
        horizontalModel.h DGL_Berechnen.c DGL_Berechnen.h Subsysteme.c Subsysteme.h Lenkmotoren.c Lenkmotoren.h matrizeCalculator.c matrizeCalculator.h)
add_executable(ansicee ${SOURCE_FILES})