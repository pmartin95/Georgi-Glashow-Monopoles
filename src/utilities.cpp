#include <iostream>
#include <algorithm>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h>
#include "sim.h"
#include "stopwatch.h"

void simulation::printAcceptance() const
{
        if(nAccepts + nRejects > 0)
        {
                std::cout << static_cast<float>(nAccepts)/ static_cast<float>(nAccepts+nRejects) * 100 << "%\n";
                std::cout << nAccepts << " accepts, " << nRejects << " rejects.\n";
        }
        else
                std::cout << "There have been no acceptances or rejections.\n";
}

void simulation::printSite(long unsigned int site_index) const
{
        for(int i =0; i < 4; i++)
                std::cout << L.site[site_index].link[i] <<std::endl;
        std::cout << L.site[site_index].higgs <<std::endl <<std::endl;
}

void simulation::printDerivatives(long unsigned int site_index) const
{
        std::cout << "Calling print derivative routine\n\n.";
        int dir;
        FORALLDIR(dir)
        std::cout << georgiGlashowActionLinkDerivative(site_index,dir, L)<< std::endl;
        std::cout << georgiGlashowActionPhiDerivative(site_index, L) << std::endl << std::endl;
        std::cout << "Exiting print derivative routine\n\n.";
}

matrix_complex CayleyHamiltonExp(const matrix_complex & A)
{

        complex<double> M = std::sqrt( A(0,1) *A(1,0) - A(0,0) *A(1,1)   );
        matrix_complex iden;
        iden.setIdentity();

        return std::cosh(M) * iden + std::sinh(M)/M *A;
}

bool isSU2(const matrix_complex & A)
{
        matrix_complex B;
        double temp;
        B = A * A.adjoint();
        temp = 2.0d - B.trace().real();
        if(std::abs(temp)>CLOSETOZERO)
                return false;
        else
                return true;
}
bool isTraceless(const matrix_complex & A)
{
        double temp;
        temp = A.trace().real();
        if(std::abs(temp)>CLOSETOZERO)
                return false;
        else
                return true;
}

bool isHermitian(const matrix_complex &A)
{
        matrix_complex B;
        double temp;
        B = A - A.adjoint();
        temp = std::abs(B(0,0)) + std::abs(B(1,0)) + std::abs(B(0,1)) + std::abs(B(1,1));
        if(std::abs(temp)>CLOSETOZERO)
                return false;
        else
                return true;
}

bool isLatticeConsistent(const lattice &L_in )
{
        long unsigned lattice_index;
        int dir_index;
        matrix_complex temp_higgs, temp_gauge[4];


        for(lattice_index = 0; lattice_index < L_in.nsites; lattice_index++)
        {
                temp_higgs = L_in.site[lattice_index].higgs;
                FORALLDIR(dir_index)
                {
                        temp_gauge[dir_index]=L_in.site[lattice_index].link[dir_index];
                        if( !isSU2(temp_gauge[dir_index]) )
                        {
                                std::cout << "Inconsistent lattice: gauge field at " << lattice_index << std::endl;
                                std::cout << temp_gauge[dir_index] << std::endl;
                        }

                }
                if( !isTraceless(temp_higgs) || !isHermitian(temp_higgs) )
                {
                        std::cout << "Inconsistent lattice: Higgs field at " << lattice_index << std::endl;
                        return false;
                }
        }

        return true;

}
