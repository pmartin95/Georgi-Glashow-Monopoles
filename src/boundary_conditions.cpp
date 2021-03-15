#include <iostream>
#include <algorithm>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h>
#include "sim.h"
#include "stopwatch.h"

void simulation::setupBoundaryConditions()
{
        char c;
        std::cout << "Enter a character to notate which type of boundary condition "
                  << "you would like to use:\n";
        std::cin >> c;
        switch (c) {
        case 'p': case 'P':
                boundary_condition = &simulation::periodicBoundaryCondition;
                break;
        default:
                std::cout << "ERROR: no valid input. Using periodic boundary conditions.\n";
                boundary_condition = &simulation::periodicBoundaryCondition;
        }
}
void simulation::setupBoundaryConditions( char boundaryType)
{
        switch (boundaryType) {
        case 'p': case 'P':
                boundary_condition = &simulation::periodicBoundaryCondition;
                break;
        default:
                std::cout << "ERROR: no valid input. Using periodic boundary conditions.\n";
                boundary_condition = &simulation::periodicBoundaryCondition;
        }
}

const matrix_complex simulation::periodicBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index,const int jump[4]) const
{
        int i;
        bool shouldBreak;
        FORALLDIR(i)
        {
                shouldBreak = (jump[i] != 0);
                if(shouldBreak)
                        break;
        }
        if(!shouldBreak)
        {
                if(matrix_num<4)
                        return L_in.site[index].link[matrix_num];
                else
                        return L_in.site[index].higgs;
        }
        int x[4];
        long unsigned int new_index;
        L_in.indexToCoordinate(index, x);
        int dir;
        FORALLDIR(dir)
        {
                x[dir] += jump[dir];
                if(x[dir] >= L_in.ns[dir] || x[dir] < 0)
                        x[dir] = (x[dir] + L_in.ns[dir])%(L_in.ns[dir]);
        }
        new_index = L_in.coordinateToIndex(x);
        if(matrix_num<4)
                return L_in.site[new_index].link[matrix_num];
        else
                return L_in.site[new_index].higgs;
}

//const matrix_complex simulation::cBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]);
//const matrix_complex simulation::twistedBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]);
