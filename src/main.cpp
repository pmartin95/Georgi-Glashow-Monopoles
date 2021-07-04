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



int main()
{

        int iters = 100;
        double energy_diffs = 0.0;
        simulation sim1;
        for (size_t i = 0; i < 50; i++)
        {
               sim1.initializeHMC();
        }
        
        
        for (size_t j = 0; j < 100; j++)
        {
                for (int i = 0; i < iters; i++)
                {
                        energy_diffs += abs(sim1.runLeapfrogSimulation());
                }
                std::cout << energy_diffs / static_cast<double>(iters) << std::endl;
                energy_diffs = 0.0;
        }

        // simulation sim1;
        // sim1.inputScheduleParameters("annealingSchedule.txt");
        // sim1.runHMCSimulationSchedule(30,300,10);
        // sim1.printDataFile("HMCtestfiles/modHMCTest1.txt");
        // sim1.stringTension("HMCtestfiles/stringTensionModHMCTest1.txt",5);

        std::cout << '\a';
        return 0;
}
