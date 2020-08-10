#pragma once
#include "lattice.h"

#ifndef _RAND_
#define _RAND_

unsigned seedGen();

const matrix_complex normalHermitianMatrix( std::mt19937_64& g ); //Also traceless
const matrix_complex uniformSU2Matrix( std::mt19937_64& g );
const double uniformReal(std::mt19937_64& g);
const double uniformReal(std::mt19937_64& g, double low, double high);

#endif
