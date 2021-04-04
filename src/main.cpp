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

        std::vector<double> testvector;
        std::mt19937_64 g(seedGen());
        std::normal_distribution<double> distribution(0.0,1.0);
        double jack_error, jack_ave;
        vectorFunc f = cosave;
        for(int i =0; i < 1000; i++)
        {
                testvector.push_back( distribution(g) + M_PI /3.0  ); // +
        }
        std::cout << "Here.\n";
        if(!computeJackknifeStatistics(testvector,f,10,jack_ave,jack_error))
                std::cout << "Jackknife succeeded.\n"
                          << "Average: " << jack_ave << std::endl
                          <<"Error: " << jack_error << std::endl;


        std::cout << '\a';
        return 0;
}

double cosave(const std::vector<double>& in)
{
        double temp;
        temp = average(in);
        return cos(temp);
}
