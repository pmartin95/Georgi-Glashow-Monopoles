#include <iostream>
#include <algorithm>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h>
#include "sim.h"
#include "stopwatch.h"

void simulation::massDiffCumSum(std::vector<double> inputData,std::vector<double> inputDataError,std::vector<double> outputData,std::vector<double> outputDataError) const
{
        std::vector<double> cumSumData, cumSumError;
        int i;

        cumSumError.reserve(inputDataError.size());
        cumSumData.reserve(inputData.size());

        cumSumData.push_back(inputData[0]);
        cumSumError.push_back(inputDataError[0]);

        for(i=0; i<inputData.size()-1; i++)
                cumSumData[i+1].push_back(cumSumData[i]+inputData[i+1]);

        outputData = cumSumData;
        for(i=0; i<inputDataError.size()-1; i++)
                cumSumError[i+1].push_back(sqrt(cumSumError[i]*cumSumError[i]+ inputDataError[i+1]  * inputDataError[i+1]  ));
        outputDataError = cumSumError;

}
