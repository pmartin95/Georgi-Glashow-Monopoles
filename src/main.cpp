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

        simulation sim1, sim2;
        double current_m2,next_m2;
        int current_data_index1;
        std::vector<double> dataCollect1, dataCollect2;

        sim1.setupBoundaryConditions('t');
        sim2.setupBoundaryConditions('c');

        sim1.inputScheduleParameters("monopole_schedule.txt");
        sim2.inputScheduleParameters("monopole_schedule.txt");

        sim1.runMHMCSimulationSchedule(20000, 100);
        sim2.runMHMCSimulationSchedule(20000, 100);
        //This collects the phi^2 data
        // for(m2s)
        // {
        //         sim1.convertDataForJackknifeMonopoleMass(current_m2,dataCollect1);
        //         sim1.convertDataForJackknifeMonopoleMass(current_m2,dataCollect2);
        //
        // }

        //Convert phi^2 data into (m2-m1) * phi^2 for jackknifing

        std::cout << '\a';
        return 0;
}
