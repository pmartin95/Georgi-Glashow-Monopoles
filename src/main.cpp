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

        std::vector<std::vector<double> > v;
        std::vector<double> b(10,0.0d);
        v.push_back(b);
        v.push_back(b);
        v[0].push_back(5.0 );
        std::cout << v[0].size() << std::endl;
        // matrix_complex iden,zero;
        // iden.setIdentity();
        // zero.setZero();
        // lattice_site site1(iden,zero);
        // lattice L_in(site1);
        // simulation sim1(L_in);
        // simulation sim1;
        // for(int i = 0; i < 100; i++)
        //         sim1.initializeHMC();
        // std::cout << sim1.CreutzRatio(1,1) << std::endl;
        // std::cout << sim1.CreutzRatio(2,2) << std::endl;


        std::cout << '\a';
        return 0;
}

double cosave(const std::vector<double>& in)
{
        double temp;
        temp = average(in);
        return cos(temp);
}
