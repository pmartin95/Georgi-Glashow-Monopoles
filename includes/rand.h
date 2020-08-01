#include <chrono>
#include <random>
#include <iostream>
#include <complex>
#include <cmath>
#include "lattice.h"

#ifndef _RAND_
#define _RAND_

unsigned seedGen();

const matrix_complex normalHermitianMatrix( URNG& g ); //Also traceless
const matrix_complex uniformSU2Matrix( URNG& g );
const double uniformReal(URNG& g);
const double uniformReal(URNG& g, double low, double high);


#endif
