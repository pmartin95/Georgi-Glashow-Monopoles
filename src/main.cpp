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



        // simulation sim1;
        // std::vector<std::vector<double> > temp;
        // double temp_ave, temp_error;
        //
        //
        //
        // sim1.inputScheduleParameters("annealingSchedule.txt");
        //
        // sim1.runMHMCSimulationSchedule(2000, 100);
        // for(int i =0; i < sim1.schedule.size(); i++)
        // {
        //         sim1.convertDataForJackknifeCreutz(sim1.schedule[i].g_value,2, temp);
        //         computeJackknifeStatistics(temp[2], average,  2, temp_ave, temp_error );
        //         std::cout << 4.0/sim1.schedule[i].g_value*sim1.schedule[i].g_value << " " << temp_ave<< " " << temp_error << std::endl;
        // }
        // sim1.resetDataPoints();
        // sim1.runHMCSimulationSchedule(30,200,10);
        // for(int i =0; i < sim1.schedule.size(); i++)
        // {
        //         sim1.convertDataForJackknifeCreutz(sim1.schedule[i].g_value,2, temp);
        //         computeJackknifeStatistics(temp[2], average,  2, temp_ave, temp_error );
        //         std::cout << 4.0/sim1.schedule[i].g_value*sim1.schedule[i].g_value << " " << temp_ave<< " " << temp_error << std::endl;
        // }

        // sim1.stringTension("string2go",10);


        std::cout << '\a';
        return 0;
}

double cosave(const std::vector<double>& in)
{
        double temp;
        temp = average(in);
        return cos(temp);
}
