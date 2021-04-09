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
        setupBoundaryConditions(c);
}

void simulation::setupBoundaryConditions( char boundaryType)
{
        switch (boundaryType) {
        case 'p': case 'P':
                boundary_condition = &simulation::periodicBoundaryCondition;
                break;
        case 'c': case 'C':
                boundary_condition = &simulation::cBoundaryCondition;
                break;
        case 't': case 'T':
                boundary_condition = &simulation::twistedBoundaryCondition;
                break;
        default:
                std::cout << "ERROR: no valid input. Using periodic boundary conditions.\n";
                boundary_condition = &simulation::periodicBoundaryCondition;
        }
}


const matrix_complex simulation::periodicBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index,const int jump[4]) const
{
        long unsigned int new_index;

        if(!isJump(jump))
                return L_in.directMatCall(index,matrix_num);
        new_index = L_in.shiftToLattice(index,jump);
        return L_in.directMatCall(new_index,matrix_num);
}



const matrix_complex simulation::cBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]) const
{
        long unsigned int new_index;
        int i,dir,x[4],y[4];
        matrix_complex temp_mat;

        if(!isJump(jump))
                return L_in.directMatCall(index,matrix_num);

        L_in.indexToCoordinate(index, x);
        FORALLDIR(dir)
        {
                x[dir] += jump[dir];
                y[dir] = L_in.shiftToLattice(x[dir],dir);
        }

        temp_mat = L_in.directMatCall(L_in.coordinateToIndex(y),matrix_num);

        for(dir = 1; dir < 4; dir++)
        {
                while(x[dir] != y[dir])
                {
                        x[dir] = L_in.incrementCoordinate(x[dir],dir);
                        temp_mat = cTwist(temp_mat,matrix_num);
                }
        }
        return temp_mat;
}

const matrix_complex simulation::cTwist(const matrix_complex& mat_A, int matrix_num) const
{
        if(matrix_num<4)
                return cTwistGaugeField(mat_A);
        else
                return cTwistHiggsField(mat_A);
}

const matrix_complex simulation::cTwistGaugeField(const matrix_complex& mat_A) const
{
        matrix_complex temp;
        temp(0,0) = mat_A(1,1);
        temp(1,1) = mat_A(0,0);
        temp(0,1) = -mat_A(1,0);
        temp(1,0) = -mat_A(0,1);
        return temp;
}

const matrix_complex simulation::cTwistHiggsField(const matrix_complex& mat_A) const
{
        return -cTwistGaugeField(mat_A);
}

const matrix_complex simulation::twistedBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]) const
{
        long unsigned int new_index;
        int i,dir,x[4],y[4];
        matrix_complex temp_mat;


        if(!isJump(jump))
                return L_in.directMatCall(index,matrix_num);

        //The following block breaks apart the index into coordinates. It then
        // creates another set of coordinates, y and shifts them onto the lattice.
        //It then find the matrix at the y coordinates and saves it to temp_mat.
        //Temp_mat is then ctwisted while x is being incremeneted to y. Note that
        //I do not increment the time direction because we are treating time periodically

        L_in.indexToCoordinate(index, x);
        FORALLDIR(dir)
        {
                x[dir] += jump[dir];
                y[dir] = L_in.shiftToLattice(x[dir],dir);
        }

        temp_mat = L_in.directMatCall(L_in.coordinateToIndex(y),matrix_num);

        for(dir = 1; dir < 4; dir++)
        {
                while(x[dir] != y[dir])
                {
                        x[dir] = L_in.incrementCoordinate(x[dir],dir);
                        temp_mat = TwistField(temp_mat,matrix_num,dir);
                }
        }
        return temp_mat;
}
const matrix_complex simulation::TwistField(const matrix_complex& mat_A,int matrix_number,int dir) const
{
        switch (dir) {
        case 0:
                return mat_A;
        case 1:
                return xTwist(mat_A,matrix_number);
        case 2:
                return yTwist(mat_A,matrix_number);
        case 3:
                return zTwist(mat_A,matrix_number);
        default:
                std::cout << "ERROR: incorrect direciton input.\n";
                exit(1);
        }
}
const matrix_complex simulation::xTwist(const matrix_complex& mat_A,int matrix_number) const
{
        matrix_complex temp;
        temp(0,0) = mat_A(1,1);
        temp(1,1) = mat_A(0,0);
        temp(0,1) = mat_A(1,0);
        temp(1,0) = mat_A(0,1);
        if(matrix_number < 4)
                return temp;
        else
                return -temp;
}
const matrix_complex simulation::yTwist(const matrix_complex& mat_A,int matrix_number) const
{
        matrix_complex temp;
        temp(0,0) = mat_A(1,1);
        temp(1,1) = mat_A(0,0);
        temp(0,1) = -mat_A(1,0);
        temp(1,0) = -mat_A(0,1);
        if(matrix_number < 4)
                return temp;
        else
                return -temp;
}
const matrix_complex simulation::zTwist(const matrix_complex& mat_A,int matrix_number) const
{
        matrix_complex temp;
        temp(0,1) = -mat_A(0,1);
        temp(1,0) = -mat_A(1,0);
        if(matrix_number < 4)
                return temp;
        else
                return -temp;
}
