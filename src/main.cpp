#include <cstdlib>
#include <sys/time.h>
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
        std::vector<unsigned int> N_Steps;
        std::vector<double> H_diff_average;
        simulation sim1;
        for(int i=0; i<30; i++)
                sim1.initializeHMC();
        int temp_step = 2;
        double temp_hdiff = 0.0d;
        int i,j;

        for(i=0; i<10; i++)
        {
                N_Steps.push_back(temp_step);
                sim1.setupSteps(temp_step);
                std::cout << "Performing " << j << " runs with " << temp_step << " intervals\n";
                for(j=0; j<30; j++)
                        temp_hdiff += abs( sim1.runLeapfrogSimulation() );
                std::cout << "step size: " << sim1.stepSize << std::endl;
                std::cout << "On to the next one!" << std::endl;
                H_diff_average.push_back(temp_hdiff/static_cast<double>(j));
                temp_step *= 2;

        }


        std::ofstream datafile;
        datafile.open("Hdiffdata.txt");
        for(i=0; i<10; i++)
                datafile << N_Steps[i] << " " << H_diff_average[i] << std::endl;
        datafile.close();
        return 0;
}
