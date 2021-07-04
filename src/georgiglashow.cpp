#include <iostream>
#include <algorithm>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h>
#include "sim.h"
#include "stopwatch.h"



////NOTE FOR PAUL
// I think the issue lies in the phiDerivativePart
//Action functions
double simulation::georgiGlashowLagrangianDensity(long unsigned int site_index) const
{
        return georgiGlashowLagrangianDensity(site_index,L);
}

double simulation::georgiGlashowLagrangianDensity(long unsigned int site_index, const lattice& L_in) const
{
        //Declaration
        matrix_complex mainPhi, mainLink[4];
        double phiSquareTrace;
        double phiDerivativePart = 0.0, plaquettePart = 0.0, miscPhiPart = 0.0;
        int jumpNone[4] = {0};
        int dir,i,j;
        //Initialization
        mainPhi = matCall(L_in,5,site_index,jumpNone);
        FORALLDIR(dir)
        mainLink[dir] = matCall(L_in,dir,site_index,jumpNone);
        phiSquareTrace = (mainPhi * mainPhi).trace().real();

        //Phi "derivative"
        phiDerivativePart = 4.0 * phiSquareTrace;
        FORALLDIR(dir)
        {
                int jumpTemp[4] = {0};
                jumpTemp[dir]++;
                phiDerivativePart -= (mainPhi * mainLink[dir] * matCall(L_in,dir,site_index,jumpTemp) * mainLink[dir].adjoint() ).trace().real();
        }
        phiDerivativePart *= 2.0 * dsv; //dsv let's me turn off this phi derivative part on and off

        //plaquette
        FORALLDIRLESSTHAN(i,j)
        {
                int muJump[4] = {0}, nuJump[4] = {0};
                muJump[i]++; nuJump[j]++;
                plaquettePart += 2.0 - ( mainLink[i] * matCall(L_in,i,site_index,muJump) * matCall(L_in,j,site_index,nuJump).adjoint() * mainLink[j].adjoint() ).trace().real();
        }
        plaquettePart *= 2.0/(g*g);
        //other phi terms
        miscPhiPart = m2 * phiSquareTrace + lambda * phiSquareTrace *phiSquareTrace;


        return phiDerivativePart + plaquettePart + miscPhiPart;
}

double simulation::georgiGlashowAction() const
{
        return georgiGlashowAction(L);
}

double simulation::georgiGlashowAction(const lattice& L_in) const
{
        double total = 0.0;
  #pragma omp parallel default(shared)
        {
                double local_total = 0.0;
    #pragma omp for
                for(long unsigned int i = 0; i<L_in.nsites; i++)
                {
                        local_total += georgiGlashowLagrangianDensity(i,L_in);
                }
    #pragma omp atomic
                total+=local_total;
        }
        return total;
}

double simulation::kineticTerm(const Plattice& P_in) const
{
        double momenta_total = 0.0;
        matrix_complex momenta_matrix;

        long unsigned site_index;
        momenta_matrix.setZero();
  #pragma omp parallel default(shared)
        {
                matrix_complex momenta_matrix_higgs,momenta_matrix_gauge;
                momenta_matrix_gauge.setZero();
                momenta_matrix_higgs.setZero();
                int i;
    #pragma omp for
                for(site_index=0; site_index < P_in.nsites; site_index++)
                {
                        FORALLDIR(i)
                        momenta_matrix_gauge+= P_in.site[site_index].link[i] *  P_in.site[site_index].link[i];
                        momenta_matrix_higgs+=   P_in.site[site_index].higgs *  P_in.site[site_index].higgs;
                        #ifdef __CHECK_NAN__
                        if( isnan(momenta_matrix_gauge.trace().real()) ||  isnan(momenta_matrix_higgs.trace().real())  )
                        {
                                std::cout << "NAN ERROR: at site " << site_index << std::endl;
                                FORALLDIR(i)
                                std::cout <<   P_in.site[site_index].link[i] << std::endl;
                                std::cout <<   P_in.site[site_index].higgs << std::endl;
                                exit(1);
                        }
                        #endif
                }

    #pragma omp critical
                {
                        momenta_matrix += momenta_matrix_gauge;
                        momenta_matrix += momenta_matrix_higgs;
                }
        }


        momenta_total = momenta_matrix.trace().real();
        return momenta_total;
}


double simulation::georgiGlashowHamiltonian(const lattice& L_in, const Plattice& P_in) const
{
        matrix_complex momenta_matrix;
        double field_total, kinetic_total,total;
        long unsigned int site_index;
        int i;

        field_total = georgiGlashowAction(L_in);
        //std::cout << "Field total: " << field_total << std::endl;

        kinetic_total = kineticTerm(P_in);
        //std::cout << "Momenta total: " << kinetic_total << std::endl;
        total = kinetic_total + field_total;
        return total;
}

double simulation::Hamiltonian() const
{
        return georgiGlashowHamiltonian(L,P) /L.nsites;
}

const matrix_complex simulation::plaquette(long unsigned site_index, int dir1, int dir2) const
{
        return plaquette(L,site_index,dir1,dir2);
}

const matrix_complex simulation::plaquette(const lattice& L_in,long unsigned site_index, int dir1, int dir2) const
{
        matrix_complex u1, u2, u3, u4;
        int jumpNone[4] ={0}, jump2[4]={0}, jump3[4]={0};
        jump2[dir1] += 1; jump3[dir2] += 1;
        u1 = matCall(L_in,dir1, site_index,jumpNone);
        u2 = matCall(L_in,dir2, site_index,jump2);
        u3 = matCall(L_in,dir1, site_index,jump3).adjoint();
        u4 = matCall(L_in,dir2, site_index,jumpNone).adjoint();

        return u1*u2*u3*u4;
}
//This function computes the mixed term in the georgi glashow lagrangian
const matrix_complex simulation::mixedGaugeHiggsTerm(const lattice& L_in, long unsigned site_index, int dir) const
{
        matrix_complex u1, u2, u3, u4;
        int jump[4] ={0},jumpNone[4] = {0};
        jump[dir]++;
        u1 = matCall(L_in,4, site_index,jumpNone);
        u2 = matCall(L_in,dir, site_index,jumpNone);
        u3 = matCall(L_in,4, site_index,jump);
        u4 = matCall(L_in,dir, site_index,jumpNone).adjoint();

        return u1*u2*u3*u4;
}
