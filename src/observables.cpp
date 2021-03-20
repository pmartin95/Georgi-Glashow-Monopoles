#include <iostream>
#include <algorithm>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h>
#include "sim.h"
#include "stopwatch.h"

//Observables
double simulation::averagePlaquettes() const
{
        int dir1,dir2;
        double subtotal = 0.0;
        unsigned long int site_index;
  #pragma omp parallel private(dir1,dir2) default(shared)
        {
                double local_subtotal = 0.0;
    #pragma omp for
                for(long unsigned int i = 0; i< L.nsites; i++)
                {
                        FORALLDIRLESSTHAN(dir1,dir2)
                        local_subtotal+= plaquette(site_index,dir1,dir2).trace().real();
                }
    #pragma omp critical
                subtotal+=local_subtotal;
        }
        return subtotal / static_cast<double>(L.nsites);
}

const matrix_complex simulation::averagePhi() const
{
        matrix_complex subtotal;
        subtotal.setZero();
        unsigned long int site_index;
  #pragma omp parallel default(shared)
        {
                matrix_complex local_subtotal;
                subtotal.setZero();
    #pragma omp for
                for(long unsigned int i = 0; i< L.nsites; i++)
                {
                        local_subtotal+= L.site[i].higgs;
                }
    #pragma omp critical
                subtotal+=local_subtotal;
        }
        return subtotal / static_cast<double>(L.nsites);
}

const matrix_complex simulation::averagePhi2() const
{
        matrix_complex subtotal;
        subtotal.setZero();
        unsigned long int site_index;
  #pragma omp parallel default(shared)
        {
                matrix_complex local_subtotal;
                subtotal.setZero();
    #pragma omp for
                for(long unsigned int i = 0; i< L.nsites; i++)
                {
                        local_subtotal += L.site[i].higgs * L.site[i].higgs;
                }
    #pragma omp critical
                subtotal += local_subtotal;
        }
        return subtotal / static_cast<double>(L.nsites);
}
// Compute the X(1,1) Creutz ratio
double simulation::CreutzRatio() const
{
        double W11;
        W11 = averagePlaquettes();
        if(W11 > 0.0)
                return -log(W11);
        else
                return 0;

}
