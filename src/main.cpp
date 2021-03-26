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
#include "stopwatch.h"
#include <unsupported/Eigen/MatrixFunctions>

int main()
{

        std::string filename = "annealingSchedule.txt";
        simulation sim1;
        sim1.switchDSV();
        sim1.inputScheduleParameters(filename);
        void runHMCSimulationSchedule(int init_iter,int iter, int iter_measure);
        sim1.runHMCSimulationSchedule(30,600,20);
        sim1.printDataFile("HMCsim1.txt");
        // sim1.resetDataPoints();
        // std::cout << "Finished running hybrid Monte Carlo.\n";
        // sim1.runMHMCSimulationSchedule(2000,100);
        // sim1.printDataFile("MHMCsim1.txt");
        return 0;
}
