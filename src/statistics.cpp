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
        return computeJackknifeStatistics(inputData,&average,setLength,Jackknife_ave,Jackknife_error);
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
