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
        //g schedule

        //std::vector<double> g_schedule{1.15,1.3,1.45,1.6,1.75,1.9,2,2.15,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3};
        std::vector<double> g_schedule;
        for(double z = 0.05; z<.75; z+=0.05)
                g_schedule.push_back(1.0/sqrt(z));
        //Hybrid Monte Carlo Simulation
        simulation sim1;
        sim1.switchDSV();
        vector<double> sim1_average_plaquettes;
        for(double g_in : g_schedule)
        {
                vector<double> sim1_plaquettes;
                sim1.setupParams(sim1.m2,sim1.lambda, g_in);
                for(int i=0; i<30; i++)
                        sim1.initializeHMC();
                for (int i = 0; i < 30; i++)
                {
                        sim1.runLeapfrogSimulation();
                        sim1_plaquettes.push_back(sim1.averagePlaquettes());
                }
                sim1_average_plaquettes.push_back(averageDoubleVector(sim1_plaquettes));
        }
        //Metropolis Hastings Simulation
        simulation sim2;
        sim2.switchDSV();
        vector<double> sim2_average_plaquettes;
        for(double g_in : g_schedule)
        {
                vector<double> sim2_plaquettes;
                sim2.setupParams(sim2.m2,sim2.lambda, g_in);
                for (int i = 0; i < 30; i++)
                {
                        sim2.multiSweepMHMC(10);
                        sim2_plaquettes.push_back(sim2.averagePlaquettes());
                }
                sim2_average_plaquettes.push_back(averageDoubleVector(sim2_plaquettes));
        }
        ;
        std::cout << "ln(g^2)   HMC    MHMC\n";
        for(int i=0; i<sim1_average_plaquettes.size(); i++)
        {
                double x,y, lng2;
                lng2 = log( g_schedule[i]   );
                x= -log(sim1_average_plaquettes[i]);
                y= -log(sim2_average_plaquettes[i]);
                std::cout << lng2 << " " << x << " " << y << std::endl;

        }

        return 0;
}
