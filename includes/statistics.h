#pragma once

#ifndef __STATISTICS__

typedef double (*vectorFunc) (const std::vector<double>& );
typedef double (*vectorVectorFunc) (const std::vector<std::vector<double>>& );

#include <cstdlib>
#include <iomanip>
#include <sys/time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <omp.h>
#include "lattice.h"
#include "rand.h"
#include "sim.h"
#include "statistics.h"
#include "stopwatch.h"
#include <unsupported/Eigen/MatrixFunctions>

#define __STATISTICS__


enum jackknife{JACKKNIFE_SUCCESS, JACKKNIFE_NO_DIVIDE};
double average(const std::vector<double>& input);
double sumVec(const std::vector<double>& input);

int computeJackknifeStatistics(const std::vector<double>& inputData, int setLength, double& Jackknife_ave, double& Jackknife_error );
int computeJackknifeStatistics(const std::vector<double>& inputData, vectorFunc f,  int setLength, double& Jackknife_ave, double& Jackknife_error );
int computeJackknifeStatistics(const std::vector<std::vector<double>>& inputData, vectorVectorFunc f,  int setLength, double& Jackknife_ave, double& Jackknife_error );
double cosave(const std::vector<double>& in);

double CreutzRatio(const std::vector<std::vector<double>>& rectangleData); 
#endif
