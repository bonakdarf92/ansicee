cmake_minimum_required(VERSION 3.5)
project(ansicee)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -lgsl -lgslcblas ")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")

#include_directories(-I/usr/include/gsl)
include_directories(include)
include_directories(src)
include_directories(Benchmark)

set(SOURCE_FILES main.c include/linear.h src/linear.c src/horizontalModel.c include/horizontalModel.h Benchmark/InitTest.c Benchmark/InitTest.h src/ReglerIO.c include/ReglerIO.h Benchmark/ErrorCorrection.c Benchmark/ErrorCorrection.h src/calcMatrix_A_B.c include/calcMatrix_A_B.h)
#set(SOURCE_FILES src/modelAntrieb.c include/modelAntrieb.h src/calcMatrix_A_B.c include/calcMatrix_A_B.h src/DGL_Berechnen.c include/DGL_Berechnen.h src/Subsysteme.c include/Subsysteme.h src/Lenkmotoren.c include/Lenkmotoren.h src/matrizeCalculator.c include/matrizeCalculator.h )
#set(SOURCE_FILES RTMakerDebug/RT_Main.c RTMakerDebug/App.c RTMakerDebug/App.h RTMakerDebug/IO.c RTMakerDebug/IO.h RTMakerDebug/IO_ADMA_User.c RTMakerDebug/IO_ADMA_User.h RTMakerDebug/IO_SLS400.c RTMakerDebug/IO_SLS400.h RTMakerDebug/IO_SLS400_User.c RTMakerDebug/IO_SLS400_User.h RTMakerDebug/app_tmp.c
        #RTMakerDebug/RTByteSwap.h RTMakerDebug/)

add_executable(ansicee ${SOURCE_FILES})