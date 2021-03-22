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

const matrix_complex simulation::directMatCall(const lattice& L_in,unsigned long int site_index,int matrix_num) const
{
        if(matrix_num<4)
                return L_in.site[site_index].link[matrix_num];
        else
                return L_in.site[site_index].higgs;
}
//The following function shifts a coordinate to be on lattice via adding and
//reducing based on lattice size
int simulation::shiftToLattice(const lattice& L_in,int coordinate, int dir) const
{

        int temp_coordinate = coordinate;
        while(temp_coordinate < 0.0)
        {
                temp_coordinate += L_in.ns[dir];
        }
        return (temp_coordinate)%(L_in.ns[dir]);
}

int simulation::incrementCoordinate(const lattice& L_in,int coordinate,int dir) const
{
        if(coordinate > L_in.ns[dir] )
                return coordinate - L_in.ns[dir];
        else if(coordinate <  0)
                return coordinate + L_in.ns[dir];
        else
        {
                std::cout << "ERROR: Indices not corresponding" << std::endl;
                std::exit(1);
        }
}

const matrix_complex simulation::periodicBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index,const int jump[4]) const
{
        long unsigned int new_index;
        int i,dir,x[4];
        bool isJump;


        //The following checks to see if there is a jump at all.
        //If not, then the pure link/Higgs field should just be returned directly
        //For clarification, if isJump is true then there is some jump
        //if isJump is false, there is no jump
        FORALLDIR(i)
        {
                isJump = (jump[i] != 0);
                if(isJump)
                        break;
        }
        if(!isJump)
                return directMatCall(L_in,index,matrix_num);

        //The following decomposes the index into coordinates and then detects
        //if the called direction is over the edge of the lattice. If so, then
        // it does a modulus on the operation that returns it to the lattice.
        //It then converts the coordinates back to an index and returns the
        //corresponding field at that site.
        L_in.indexToCoordinate(index, x);
        FORALLDIR(dir)
        {
                x[dir] += jump[dir];
                if(x[dir] >= L_in.ns[dir] || x[dir] < 0)
                        x[dir] = shiftToLattice(L_in,x[dir],dir);
        }
        new_index = L_in.coordinateToIndex(x);
        return directMatCall(L_in, new_index,matrix_num);
}



const matrix_complex simulation::cBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]) const
{
        long unsigned int new_index;
        int i,dir,x[4],y[4];
        bool isJump;
        matrix_complex temp_mat;
        //The following checks to see if there is a jump at all.
        //If not, then the pure link/Higgs field should just be returned directly
        //For clarification, if isJump is true then there is some jump
        //if isJump is false, there is no jump
        FORALLDIR(i)
        {
                isJump = (jump[i] != 0);
                if(isJump)
                        break;
        }
        if(!isJump)
                return directMatCall(L_in,index,matrix_num);

        //The following block breaks apart the index into coordinates. It then
        // creates another set of coordinates, y and shifts them onto the lattice.
        //It then find the matrix at the y coordinates and saves it to temp_mat.
        //Temp_mat is then ctwisted while x is being incremeneted to y. Note that
        //I do not increment the time direction because we are treating time periodically

        L_in.indexToCoordinate(index, x);
        FORALLDIR(dir)
        {
                x[dir] += jump[dir];
                y[dir] = shiftToLattice(L_in,x[dir],dir);
        }

        temp_mat = directMatCall(L_in,L_in.coordinateToIndex(y),matrix_num);

        for(dir = 1; dir < 4; dir++)
        {
                while(x[dir] != y[dir])
                {
                        x[dir] = incrementCoordinate(L_in,x[dir],dir);
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
        bool isJump;
        matrix_complex temp_mat;
        //The following checks to see if there is a jump at all.
        //If not, then the pure link/Higgs field should just be returned directly
        //For clarification, if isJump is true then there is some jump
        //if isJump is false, there is no jump
        FORALLDIR(i)
        {
                isJump = (jump[i] != 0);
                if(isJump)
                        break;
        }
        if(!isJump)
                return directMatCall(L_in,index,matrix_num);

        //The following block breaks apart the index into coordinates. It then
        // creates another set of coordinates, y and shifts them onto the lattice.
        //It then find the matrix at the y coordinates and saves it to temp_mat.
        //Temp_mat is then ctwisted while x is being incremeneted to y. Note that
        //I do not increment the time direction because we are treating time periodically

        L_in.indexToCoordinate(index, x);
        FORALLDIR(dir)
        {
                x[dir] += jump[dir];
                y[dir] = shiftToLattice(L_in,x[dir],dir);
        }

        temp_mat = directMatCall(L_in,L_in.coordinateToIndex(y),matrix_num);

        for(dir = 1; dir < 4; dir++)
        {
                while(x[dir] != y[dir])
                {
                        x[dir] = incrementCoordinate(L_in,x[dir],dir);
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
