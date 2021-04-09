#include <cstdlib>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <omp.h>
#include "statistics.h"
#include <unsupported/Eigen/MatrixFunctions>

double average(const std::vector<double>& input)
{
        double temp_sum = 0.0d;
        int vec_size = input.size();

        for(int i=0; i < vec_size; i++)
                temp_sum += input[i];

        return temp_sum / static_cast<double>(vec_size);

}
double sumVec(const std::vector<double>& input)
{
        double temp_sum = 0.0d;
        int vec_size = input.size();

        for(int i=0; i < vec_size; i++)
                temp_sum += input[i];

        return temp_sum;
}
int computeJackknifeStatistics(const std::vector<double>& inputData, int setLength, double& Jackknife_ave, double& Jackknife_error )
{
        return computeJackknifeStatistics(inputData,average,setLength,Jackknife_ave,Jackknife_error);
}

int computeJackknifeStatistics(const std::vector<double>& inputData, vectorFunc f,  int setLength, double& Jackknife_ave, double& Jackknife_error )
{
        if(inputData.size()%setLength != 0)
                return JACKKNIFE_NO_DIVIDE;
        std::vector<double> holdingVector;

        int numSets;
        numSets = inputData.size()/setLength;
        std::vector<double> setOfAverages(numSets,0.0);
        for(int i = 0; i < numSets; i++ )
        {
                holdingVector.clear();
                holdingVector.reserve(inputData.size()-setLength);
                holdingVector.insert(holdingVector.end(),inputData.begin(),inputData.begin() +i*setLength  );
                holdingVector.insert(holdingVector.end(),inputData.begin() + (i+1)*setLength,inputData.end());
                setOfAverages[i]= f(holdingVector);
        }
        Jackknife_ave = average(setOfAverages);

        std::vector<double> deviations(numSets,0.0);
        double temp;
        for(int i = 0; i < numSets; i++)
        {
                temp = setOfAverages[i] - Jackknife_ave;
                deviations[i] = temp * temp;
        }
        Jackknife_error = sqrt( static_cast<double>(numSets-1) * average(deviations)  );
        return JACKKNIFE_SUCCESS;
}


int computeJackknifeStatistics(const std::vector<std::vector<double> >& inputData, vectorVectorFunc f,  int setLength, double& Jackknife_ave, double& Jackknife_error )
{
        if(inputData[0].size()%setLength != 0)
                return JACKKNIFE_NO_DIVIDE;
        std::vector<std::vector<double> > holdingVector;
        holdingVector.reserve(inputData.size());
        std::vector<double> v = {0.0};
        for(int i =0; i<inputData.size(); i++ )
                holdingVector.push_back(v);

        int numSets;
        numSets = inputData[0].size()/setLength;
        std::vector<double> setOfAverages(numSets,0.0);
        for(int i = 0; i < numSets; i++ )
        {
                for(int j = 0; j < inputData.size(); j++)
                {
                        holdingVector[j].erase(holdingVector[j].begin(), holdingVector[j].end());
                        holdingVector[j].reserve(inputData[j].size()-setLength);
                        holdingVector[j].insert(holdingVector[j].end(),inputData[j].begin(),inputData[j].begin() +i*setLength  );
                        holdingVector[j].insert(holdingVector[j].end(),inputData[j].begin() + (i+1)*setLength,inputData[j].end());
                }

                setOfAverages[i]= f(holdingVector);
        }
        Jackknife_ave = average(setOfAverages);

        std::vector<double> deviations(numSets,0.0);
        double temp;
        for(int i = 0; i < numSets; i++)
        {
                temp = setOfAverages[i] - Jackknife_ave;
                deviations[i] = temp * temp;
        }
        Jackknife_error = sqrt( static_cast<double>(numSets-1) * average(deviations)  );
        return JACKKNIFE_SUCCESS;
}



//First vector is W(R,R), then W(R,R-1), then W(R-1,R-1)
double CreutzRatio(const std::vector<std::vector<double> >& rectangleData)
{
        if( rectangleData.size() != 3)
        {
                std::cout << "Data error: improper vector size.\n"
                          << "Actual size is " << rectangleData.size() << std::endl;
                exit(1);

        }

        double WRR, WRR1, WR1R1,temp;
        WRR = average(rectangleData[0]);
        WRR1 = average(rectangleData[1]);
        WR1R1 = average(rectangleData[2]);
        WRR1 = WRR1 * WRR1;
        temp = WRR * WR1R1/WRR1;
        if(temp>0.0)
                return -log(temp);
        else
                return 0.0;
}
