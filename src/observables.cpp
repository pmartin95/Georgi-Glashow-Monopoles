#include <iostream>
#include <algorithm>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h>
#include "sim.h"
#include "stopwatch.h"

double simulation::averageWilsonRectangle(int dir1_len,int dir2_len) const
{
        double cumulative_value = 0.0;
        int i, j;
        for(unsigned long site_index=0; site_index < L.nsites; site_index++) FORALLDIR(i) FORALLDIRBUT(j,i)
                {
                        cumulative_value += rectangleWilson(site_index,i,dir1_len,j,dir2_len);
                }
        return cumulative_value/static_cast<double>(4*3*L.nsites);
}

double simulation::rectangleWilson(unsigned long site_index, int dir1,int dir1_len, int dir2, int dir2_len) const
{

        matrix_complex cumulative;
        long unsigned current_site_index;
        int current_coordinates[4] = {0};
        int jump_coordinates[4] = {0};
        int i;

        cumulative.setIdentity();
        current_site_index = site_index;
        L.indexToCoordinate(site_index,current_coordinates);
        if(dir1_len == 0 || dir2_len == 0)
                return 1.0d;
        //Bottom part of Wilson rectangle
        for(i=0; i < dir1_len; i++)
        {
                cumulative *= matCall(L,dir1,site_index,jump_coordinates);
                current_coordinates[dir1]++;
                jump_coordinates[dir1]++;
                current_site_index = L.coordinateToIndex(current_coordinates);
        }
        //Right part of Wilson rectangle
        for(i=0; i < dir2_len; i++)
        {

                cumulative *= matCall(L,dir2,site_index,jump_coordinates);
                current_coordinates[dir2]++;
                jump_coordinates[dir2]++;
                current_site_index = L.coordinateToIndex(current_coordinates);
        }
        current_coordinates[dir1]--;
        jump_coordinates[dir1]--;
        current_site_index = L.coordinateToIndex(current_coordinates);
        //Top part of Wilson rectangle
        for(i=0; i < dir1_len; i++)
        {

                cumulative *= matCall(L,dir1,site_index,jump_coordinates).adjoint();
                current_coordinates[dir1]--;
                jump_coordinates[dir1]--;
                current_site_index = L.coordinateToIndex(current_coordinates);
        }
        current_coordinates[dir1]++;
        jump_coordinates[dir1]++;
        current_coordinates[dir2]--;
        jump_coordinates[dir2]--;
        current_site_index = L.coordinateToIndex(current_coordinates);
        //Left part of Wilson rectangle
        for(i=0; i < dir2_len; i++)
        {

                cumulative *= matCall(L,dir2,site_index,jump_coordinates).adjoint();
                current_coordinates[dir2]--;
                jump_coordinates[dir2]--;
                current_site_index = L.coordinateToIndex(current_coordinates);
        }
        current_coordinates[dir2]++;
        jump_coordinates[dir2]++;
        current_site_index = L.coordinateToIndex(current_coordinates);
        //This is a check to make sure the rectangle gets back to the start
        if(!( current_site_index ==site_index    ))
                std::cout << "index did not line up.\n";
        return cumulative.trace().real()*0.5;
}
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
        return subtotal/static_cast<double>(4*3*L.nsites);
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
        return subtotal/static_cast<double>(4*3*L.nsites);
}
// // Compute the X(1,1) Creutz ratio
// double simulation::CreutzRatio() const
// {
//         double W11;
//         W11 = averagePlaquettes();
//         if(W11 > 0.0)
//                 return -log(W11);
//         else
//                 return 0;
//
// }
