cmake_minimum_required(VERSION 3.14)
project(Georgi_Glashow_monopoles CXX)

file(GLOB SRC  "src/*.cpp" )

string(APPEND CMAKE_CXX_FLAGS "-I ./../eigen -I ./../includes -lm -g -fopenmp -O3")

add_executable(GGM ${SRC})
