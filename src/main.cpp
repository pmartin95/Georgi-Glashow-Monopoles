#include <cstdlib>
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
#include "stopwatch.h"
#include <unsupported/Eigen/MatrixFunctions>

int main()
{
        std::string filename = "testsched.txt";
        std::vector<double> g(10),m2(10),lambda(10);
        simulation sim1;
        sim1.generateScheduleFile(filename,g,m2,lambda);
        return 0;
}
