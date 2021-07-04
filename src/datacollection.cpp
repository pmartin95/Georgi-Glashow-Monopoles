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
#include "statistics.h"
#include "stopwatch.h"

void simulation::appendDataPoint()
{
        data_point_t temp;
  #if RECORD_G
        temp.g_value = g;
  #endif
  #if RECORD_M2
        temp.m2_value = m2;
  #endif
  #if RECORD_LAMBDA
        temp.lambda_value = lambda;
  #endif
  #if RECORD_PHI2
        temp.average_phi2_value = averagePhi2().trace().real();
  #endif
  #if RECORD_RECTANGLES
        for(int i =0; i < RECT_SIZE; i++)
                temp.upper_rectangles[i] = averageWilsonRectangle(i,i);
        for(int j = 0; j< RECT_SIZE-1; j++)
                temp.lower_rectangles[j] =averageWilsonRectangle(j,j-1);
  #endif

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
        datafile << "#g m^2 lambda average_phi^2 upper_rectangles lower_rectangles";
        for(int i=0; i<data.size(); i++)
        {
                datafile << std::endl;
                #if RECORD_G
                datafile << std::setw(width) << data[i].g_value;
                #endif
                #if RECORD_M2
                datafile << std::setw(width) << data[i].m2_value;
                #endif
                #if RECORD_LAMBDA
                datafile << std::setw(width) << data[i].lambda_value;
                #endif
                #if RECORD_PHI2
                datafile << std::setw(width) << data[i].average_phi2_value;
                #endif
                #if RECORD_RECTANGLES
                for(int i =0; i < RECT_SIZE; i++)
                        datafile << std::setw(width) << data[i].upper_rectangles[i];
                for(int j = 0; j< RECT_SIZE-1; j++)
                        datafile << std::setw(width) << data[i].lower_rectangles[j];
                #endif
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

void simulation::reverseSchedule()
{
        vector<schedule_element_t> temp_sched(schedule.size());
        for(int i=0; i < schedule.size(); i++)
                temp_sched[i] = schedule[schedule.size()-i-1];
        schedule = temp_sched;
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
                std::cout << "Running HMC sim " << sched_index+1 << "/" << schedule.size() <<std::endl;
                for(int init_index=0; init_index<init_iter; init_index++)
                        initializeHMC();
                std::cout << "Lattice initalized.\n";
                for (int iter_index = 0; iter_index < iter; iter_index++)
                {
                        std::cout << "Running leapfrog " << iter_index+1 << "/" << iter << std::endl
                                  << " in segment " << sched_index+1 << "/" << schedule.size() <<std::endl;
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
                std::cout << "Running MHMC sim " << sched_index+1 << "/" << schedule.size() <<std::endl;
                for(int iter_index = 0; iter_index < num_iter; iter_index++)
                {
                        multiSweepMHMC(iter_measure);
                        appendDataPoint();
                }
        }
}

//Looks for values of g_in in the dataset
//Pulls out the W(R,R), W(R-1,R), and W(R-1,R-2) Wilson rectangles
void simulation::convertDataForJackknifeCreutz(double g_in,int R, std::vector<std::vector<double> >& rectData) const
{
        std::vector<double> WR[3];
        for(int i = 0; i < data.size(); i++)
                if(g_in == data[i].g_value)
                { 
                        WR[0].push_back(data[i].upper_rectangles[R]);
                        WR[1].push_back(data[i].lower_rectangles[R-1]);
                        WR[2].push_back(data[i].upper_rectangles[R-1]);
                }
        rectData.clear();
        rectData.push_back(WR[0]);
        rectData.push_back(WR[1]);
        rectData.push_back(WR[2]);
}


void simulation::convertDataForJackknifeMonopoleMass(double m2_in, std::vector<double> & phi2data) const
{
        phi2data.clear();
        for(int i = 0; i < data.size(); i++)
                if(m2_in == data[i].m2_value)
                {
                        phi2data.push_back(data[i].average_phi2_value);
                }

}

void simulation::stringTension(const std::string& base_filename,int blk_len)
{
        //go through data and compile data together that has the same g value.
        //bin the data into a vector of vectors
        //ship the binned data through Jackknife with the Creutz ratio function
        //compile g values, string tension,and error all into vectors
        //flip around g into beta
        //output the beta, tension, and error to file
        double current_g;
        int current_data_index;
        std::vector<std::vector<double> > dataCollect[RECT_SIZE-1];
        for(int i=0; i < RECT_SIZE-1; i++)
                dataCollect[i].resize(3);
        std::vector<double> final_g, betas;
        std::vector<std::vector<double> > CreutzRatios(RECT_SIZE-1), CreutzRatioErrors(RECT_SIZE-1);

        current_data_index = 0;
        while(current_data_index < data.size())
        {
                current_g = data[current_data_index].g_value;
                final_g.push_back(current_g);
                for(int R = 1; R<RECT_SIZE; R++)
                {
                        convertDataForJackknifeCreutz(current_g,R,dataCollect[R-1]);
                        double temp_ave, temp_error;
                        if( !computeJackknifeStatistics(dataCollect[R-1], CreutzRatio,  blk_len, temp_ave, temp_error ) )
                        {

                                CreutzRatios[R-1].push_back(temp_ave);
                                CreutzRatioErrors[R-1].push_back(temp_error);
                        }

                }
                while(current_g == data[current_data_index].g_value && current_data_index < data.size())
                        current_data_index++;

        }
        for(int i= 0; i < final_g.size(); i++)
                betas.push_back(4.0/final_g[i]/final_g[i]);

        ofstream file[RECT_SIZE-1];
        for(int i = 0; i < RECT_SIZE-1; i++)
        {
                std::ofstream file;
                file.open(base_filename+to_string(i)+".txt",std::ios::out);
                std::cout << "g size " << final_g.size() << std::endl;
                file << std::fixed << std::setprecision(8);
                for(int j=0; j<final_g.size(); j++)
                        if(CreutzRatios[i][j]>0.0)
                                file << betas[j]
                                     << " " << CreutzRatios[i][j]
                                     << " " <<CreutzRatioErrors[i][j] << std::endl;
                file.close();
        }


}
