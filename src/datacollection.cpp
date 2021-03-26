#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <complex>
#include <vector>

#include <string>
#include <math.h>
#include <stdlib.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h>
#include "sim.h"
#include "stopwatch.h"

void simulation::appendDataPoint()
{
        data_point_t temp;
        temp.g_value = g;
        temp.m2_value = m2;
        temp.lambda_value = lambda;
        temp.average_plaquette_value =averagePlaquettes();
        temp.average_phi2_value = averagePhi2().trace().real();
        data.push_back(temp);
}

void simulation::resetDataPoints()
{
        data.clear();
}


void simulation::printDataFile(const std::string& filename) const
{
        int width = 20;
        std::ofstream datafile;
        datafile.open(filename,std::ios::out);
        datafile << "#g m^2 lambda average_plaquette average_phi^2\n";
        for(int i=0; i<data.size(); i++)
        {
                datafile << std::setw(width) << data[i].g_value;
                datafile << std::setw(width) << data[i].m2_value;
                datafile << std::setw(width) << data[i].lambda_value;
                datafile << std::setw(width) << data[i].average_plaquette_value;
                datafile << std::setw(width) << data[i].average_phi2_value;
                datafile << std::endl;
        }
        datafile.close();
}

void simulation::inputScheduleParameters(const std::string& filename)
{
        std::ifstream inputfile;
        schedule_element_t temp_sched;
        inputfile.open(filename);
        double temp_g, temp_m2, temp_lambda;
        while(inputfile >> temp_sched.g_value >> temp_sched.m2_value >>  temp_sched.lambda_value)
        {
                schedule.push_back(temp_sched);
        }
}


void simulation::generateScheduleFile(const std::string& filename,const std::vector<double>& gs,const std::vector<double>& m2s,const std::vector<double>& lambdas)
{
        std::ofstream datafile;
        datafile.open(filename,std::ios::out);
        if(  !( (gs.size() ==  m2s.size()) || (gs.size() ==  lambdas.size()) || (m2s.size() ==  lambdas.size())  )   )
                exit(1);

        for (int i = 0; i < gs.size(); i++)
        {
                datafile << gs[i] << " " << m2s[i] << " " << lambdas[i] << " " << std::endl;
        }

}

void simulation::loadScheduleValues(int i)
{
        setupParams(schedule[i].m2_value,schedule[i].lambda_value, schedule[i].g_value);
}

void simulation::runHMCSimulationSchedule(int init_iter,int iter, int iter_measure)
{
        for(int sched_index=0; sched_index < schedule.size(); sched_index++ )
        {
                loadScheduleValues(sched_index);
                for(int init_index=0; init_index<init_iter; init_index++)
                        initializeHMC();
                for (int iter_index = 0; iter_index < iter; iter_index++)
                {
                        runLeapfrogSimulation();
                        if(iter_index%iter_measure == 0)
                                appendDataPoint();
                }
        }
}

void simulation::runMHMCSimulationSchedule(int iter, int iter_measure)
{
        int num_iter = iter/iter_measure;
        std::cout << "num_iter= " << num_iter << std::endl;
        for(int sched_index=0; sched_index < schedule.size(); sched_index++ )
        {
                loadScheduleValues(sched_index);
                for(int iter_index = 0; iter_index < num_iter; iter_index++)
                {
                        multiSweepMHMC(iter_measure);
                        appendDataPoint();
                }
        }
}
