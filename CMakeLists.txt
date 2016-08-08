cmake_minimum_required(VERSION 3.5)
project(ansicee)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/include/gsl -lgsl")

include_directories(-I/usr/include/gsl)
include_directories(include)
include_directories(src)
include_directories(Benchmark)

set(SOURCE_FILES main.c include/linear.h src/linear.c src/modelAntrieb.c include/modelAntrieb.h src/calcMatrix_A_B.c include/calcMatrix_A_B.h src/horizontalModel.c
        include/horizontalModel.h src/DGL_Berechnen.c include/DGL_Berechnen.h src/Subsysteme.c include/Subsysteme.h src/Lenkmotoren.c include/Lenkmotoren.h src/matrizeCalculator.c include/matrizeCalculator.h Benchmark/InitTest.c Benchmark/InitTest.h)
add_executable(ansicee ${SOURCE_FILES})