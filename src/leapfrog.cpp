#include <iostream>
#include <algorithm>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h>
#include "sim.h"
#include "stopwatch.h"

void simulation::updateFields(long unsigned site_index,double time_step, const Plattice & P_in,const lattice & L_in, lattice& L_out )
{
        int dir;

        #ifdef __GAUGE_EVOLUTION__
        FORALLDIR(dir)
        L_out.site[site_index].link[dir] = CayleyHamiltonExp((complex<double>(0.0,1.0) * time_step * P_in.site[site_index].link[dir] )) * L_in.site[site_index].link[dir];
        #endif

        #ifdef __HIGGS_EVOLUTION__
        L_out.site[site_index].higgs = L_in.site[site_index].higgs + time_step * P_in.site[site_index].higgs;
        #endif
}

void simulation::updateMomenta(long unsigned site_index,double time_step, const Plattice & P_in,const lattice & L_in, Plattice& P_out )
{
        int dir;

        #ifdef __GAUGE_EVOLUTION__
        FORALLDIR(dir)
        P_out.site[site_index].link[dir] = P_in.site[site_index].link[dir] - time_step * georgiGlashowActionLinkDerivative(site_index, dir, L_in);
        #endif

        #ifdef __HIGGS_EVOLUTION__
        P_out.site[site_index].higgs = P_in.site[site_index].higgs - time_step * georgiGlashowActionPhiDerivative(site_index, L_in);
        #endif
}

void simulation::initializeHMC()
{
        leapfrogOneStep();
        L = Ltemp[0];
        P = Ptemp[0];
        Ptemp[1] = Ptemp[0];
        Ltemp[1] = Ltemp[0];
}

double simulation::runLeapfrogSimulation()
{
        double Hdiff;

        leapfrogOneStep();
        Hdiff = georgiGlashowHamiltonian(L, P) - georgiGlashowHamiltonian(Ltemp[0], Ptemp[0]);
        if(metropolisDecision())
                copyLatticeAndRefresh(Ltemp[0],L);
        else
                copyLatticeAndRefresh(L,Ltemp[0]);
        Ltemp[1] = L;
        #ifdef __CHECK_LATTICE__
        isLatticeConsistent(L);
        #endif
        return Hdiff;
}

void simulation::copyLatticeAndRefresh(const lattice & L_in,lattice & L_out)
{
        L_out = L_in;
        resetMomenta();
}

void simulation::copyLatticePlattice(const lattice & L_in, const Plattice & P_in, lattice & L_out, Plattice & P_out)
{
        L_out = L_in;
        P_out = P_in;
}

void simulation::leapfrogOneStep()
{
        long unsigned int site_index;
        int i;

//Initial Step
        #pragma omp parallel private(i)
        {

                #pragma omp for
                for(site_index=0; site_index < L.nsites; site_index++)
                {
                        updateMomenta(site_index,stepSize/2.0, Ptemp[0],Ltemp[0], Ptemp[1]);
                }

//Intermediate Steps
                for(i = 0; i<steps-1; i++)
                {
                        wholeStepEvolve(Ltemp[i%2],Ptemp[(i+1)%2],  Ltemp[(i+1)%2], Ptemp[i%2] );
                        #pragma omp barrier
                        #pragma omp single
                        {
                                if( i == steps-2)
                                        copyLatticePlattice(  Ltemp[(i+1)%2],Ptemp[i%2], Ltemp[i%2],Ptemp[(i+1)%2] );
                        }
                }
//Last step
                #pragma omp for
                for(site_index=0; site_index < L.nsites; site_index++)
                        updateFields(site_index,stepSize,Ptemp[0], Ltemp[0], Ltemp[1]);

                #pragma omp for
                for(site_index=0; site_index < L.nsites; site_index++)
                        updateMomenta(site_index,stepSize/2.0, Ptemp[0],Ltemp[1],Ptemp[1]);
        }
        copyLatticePlattice(Ltemp[1],Ptemp[1],Ltemp[0],Ptemp[0]);

}

void simulation::wholeStepEvolve(lattice L_in, Plattice P_in, lattice L_out, Plattice P_out)
{
        #pragma omp for
        for(long unsigned site_index=0; site_index < L.nsites; site_index++)
        {
                //std::cout << "thread num " << omp_get_thread_num() << " is updating site " << site_index << std::endl;
                updateFields(site_index,stepSize,P_in,L_in, L_out );
        }


        #pragma omp for
        for(long unsigned site_index=0; site_index < L.nsites; site_index++)
        {
                updateMomenta(site_index,stepSize, P_in,L_in,P_out);
        }
}
//This function should return true if the new config is to be accepted
bool simulation::metropolisDecision()
{
        double expResult;
        double randomDecider;
        double H_new, H_old;
        bool updateStatus;
        H_new = georgiGlashowHamiltonian(Ltemp[0],Ptemp[0]);
        H_old = georgiGlashowHamiltonian(L,P);
        expResult = std::min(exp(H_old - H_new), 1.0);
        randomDecider = uniformReal(randomGenerator);
        updateStatus = (expResult > randomDecider);
        AcceptanceCounter(updateStatus);
        return updateStatus;
}


const matrix_complex simulation::georgiGlashowActionLinkDerivative(long unsigned int site_index, int dir, const lattice& L_in) const //
{
        matrix_complex topStaple, bottomStaple;
        matrix_complex temp1,temp2;
        matrix_complex Ulink;
        matrix_complex subtotal1,subtotal2;
        matrix_complex identityMatrix;

        int jump1[4] = {0};
        int jumpNone[4] = {0,0,0,0};
        int temp_dir;
        identityMatrix.setIdentity();

        Ulink = matCall(L_in,dir,site_index,jumpNone);
//  std::cout << "Ulink: " << Ulink << std::endl;
        topStaple.setZero(); bottomStaple.setZero();
        temp1.setZero(); temp2.setZero();
        subtotal1.setZero(); subtotal2.setZero();
        jump1[dir]++;

        FORALLDIRBUT(temp_dir,dir)
        {
                int jump2[4] = {0}, jump3[4] = {0}, jump4[4] = {0};
                jump2[temp_dir]++; jump3[temp_dir]--;
                jump4[dir]++; jump4[temp_dir]--;

                topStaple.noalias() += matCall(L_in, site_index,temp_dir,jump1) * matCall(L_in, site_index,dir,jump2).adjoint() * matCall(L_in, site_index,temp_dir,jumpNone).adjoint();
                bottomStaple.noalias() += matCall(L_in, site_index,temp_dir,jump4).adjoint() * matCall(L_in, site_index,dir,jump3).adjoint() * matCall(L_in, site_index,temp_dir,jump3);
        }

        subtotal1.noalias() += complex<double>(0.0,1.0/(g*g)) * (  Ulink*(topStaple + bottomStaple)  - (topStaple + bottomStaple).adjoint()  * Ulink.adjoint()   );
        subtotal1 = subtotal1 - subtotal1.trace() * 0.5 * identityMatrix;
        temp1 = Ulink * matCall(L_in,4,site_index,jump1) * matCall(L_in,dir,site_index,jumpNone).adjoint();

        temp2 = matCall(L_in,4,site_index,jumpNone);

        subtotal2.noalias() += (temp1*temp2 - temp2.adjoint()*temp1.adjoint()) * complex<double>(0.0,0.5*dsv);


        return subtotal1 + subtotal2;// best so far
}

const matrix_complex simulation::georgiGlashowActionPhiDerivative(long unsigned int site_index, const lattice& L_in) const
{
        matrix_complex temp1, temp2, temp3,Iden,temp;
        int temp_dir;
        int jumpNone[4] = {0,0,0,0};
        double traced_part;
        temp2.setZero();
        temp3.setZero();
        Iden.setIdentity();

        FORALLDIR(temp_dir)
        {
                int jump1[4] = {0}, jump2[4] ={0};
                jump1[temp_dir]++; jump2[temp_dir]--;
                temp2 += matCall(L_in,temp_dir,site_index,jumpNone) * matCall(L_in,4,site_index,jump1) * matCall(L_in,temp_dir,site_index,jumpNone).adjoint();
                temp3 += matCall(L_in,temp_dir,site_index,jump2).adjoint() *  matCall(L_in,4,site_index,jump2)* matCall(L_in,temp_dir,site_index,jump2);
        }
        temp1 = matCall(L_in,4,site_index,jumpNone);
        temp1 = (2.0 * lambda *temp1*(temp1*temp1).trace() + (m2+8.0*dsv)*temp1  );

        //std::cout << "Link force, temp1: " << temp1 << std::endl;
        temp = dsv*(temp2 + temp3);
        temp =temp1 - temp + temp.trace() *0.5 * Iden;
        return temp;
}
