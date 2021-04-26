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
        #ifdef __GAUGE_EVOLUTION__
        int dir;
        FORALLDIR(dir)
        L_out.site[site_index].link[dir] = CayleyHamiltonExp((complex<double>(0.0,1.0) * time_step * P_in.site[site_index].link[dir] )) * L_in.site[site_index].link[dir];
        #endif

        #ifdef __HIGGS_EVOLUTION__
        L_out.site[site_index].higgs = L_in.site[site_index].higgs + time_step * P_in.site[site_index].higgs;
        #endif
}

void simulation::updateMomenta(long unsigned site_index,double time_step, const Plattice & P_in,const lattice & L_in, Plattice& P_out )
{
        #ifdef __GAUGE_EVOLUTION__
        int dir;
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
        resetMomenta();
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


matrix_complex simulation::georgiGlashowActionLinkDerivative(long unsigned int site_index, int dir, const lattice& L_in) const //
{

        matrix_complex temp;
        temp.noalias() = georgiGlashowActionPureGaugeDerivative(site_index,dir,L_in);
  #ifdef __HIGGS_GAUGE_MIXED_TERM__
        temp.noalias() += georgiGlashowActionMixedGaugeDerivative(site_index,dir,L_in);
  #endif

        return temp;
}

matrix_complex simulation::georgiGlashowActionPureGaugeDerivative(long unsigned int site_index, int dir, const lattice& L_in) const
{
        matrix_complex Ulink,identityMatrix,subtotal,total;

        int jump0[4] = {0};
        identityMatrix.setIdentity();
        Ulink = matCall(L_in, site_index,dir,jump0);

        subtotal.noalias() = Ulink * staple(site_index,dir,L_in);
        total.noalias() =  complex<double>(0.0,0.5*invg*invg) *(  subtotal - subtotal.adjoint());
        //Just to make sure it's really traceless. In theory, total should already be traceless
        return total - total.trace() *0.5 *identityMatrix;
}

matrix_complex simulation::staple(long unsigned int site_index, int dir, const lattice& L_in) const
{
        matrix_complex topStaple, bottomStaple;
        int temp_dir;
        int jump0[4] = {0};
        topStaple.setZero();
        bottomStaple.setZero();

        FORALLDIRBUT(temp_dir,dir)
        {
                int jump1[4] = {0};
                int jump2[4] = {0};
                int jump3[4] = {0};
                int jump4[4] = {0};

                jump1[dir]++;
                jump2[temp_dir]++;
                jump3[dir]++;
                jump3[temp_dir]--;
                jump4[temp_dir]--;

                topStaple.noalias()+=matCall(L_in, site_index,temp_dir,jump1) * matCall(L_in, site_index,dir,jump2).adjoint() * matCall(L_in, site_index,temp_dir,jump0).adjoint();
                bottomStaple.noalias()+= matCall(L_in, site_index,temp_dir,jump3).adjoint() * matCall(L_in, site_index,dir,jump4).adjoint() * matCall(L_in, site_index,temp_dir,jump4);
        }
        return topStaple + bottomStaple;
}

matrix_complex simulation::georgiGlashowActionMixedGaugeDerivative(long unsigned int site_index, int dir, const lattice& L_in) const
{
        int jump0[4] = {0};
        int jump1[4] = {0};
        jump1[dir]++;
        matrix_complex Ulink, temp1, temp2;
        matrix_complex identityMatrix;
        identityMatrix.setIdentity();

        Ulink = matCall(L_in, site_index,dir,jump0);
        temp1.noalias() = Ulink * matCall(L_in, site_index,4,jump1) * Ulink.adjoint()*matCall(L_in, site_index,4,jump0);
        temp2.noalias() = complex<double>(0.0,0.5) * (temp1 - temp1.adjoint());
        //Just to make sure it's really traceless. In theory, total should already be traceless
        return temp2 - 0.5*identityMatrix * temp2.trace();
}

matrix_complex simulation::georgiGlashowActionPhiDerivative(long unsigned int site_index, const lattice& L_in) const
{
        matrix_complex temp;
        temp.noalias() = georgiGlashowActionPhiMPart(site_index,L_in) + georgiGlashowActionLambdaPart(site_index,L_in);
  #ifdef __HIGGS_GAUGE_MIXED_TERM__
        temp.noalias() += georgiGlashowActionPhiKineticPart(site_index,L_in);
  #endif
        return temp;
}

matrix_complex simulation::georgiGlashowActionPhiKineticPart(long unsigned int site_index, const lattice& L_in) const
{
        matrix_complex first, second, U1,U2, phi,temp,Iden;
        int dir;
        int jump0[4] = {0};
        Iden.setIdentity();
        phi = matCall(L_in,4,site_index,jump0);
        FORALLDIR(dir)
        {
                int jump1[4] = {0};
                int jump2[4] = {0};
                jump1[dir]++;
                jump2[dir]--;
                U1 = matCall(L_in,dir,site_index,jump0);
                U2 = matCall(L_in,dir,site_index,jump2);
                first.noalias() += U1 * matCall(L_in,4,site_index,jump1) * U1.adjoint();
                second.noalias() += U2.adjoint() * matCall(L_in,4,site_index,jump2) * U2;

        }
        temp =  8.0 * phi - (first + second);
        return temp - temp.trace() * 0.5 * Iden;
}

matrix_complex simulation::georgiGlashowActionPhiMPart(long unsigned int site_index, const lattice& L_in) const
{
        int jump0[4] = {0};
        return m2 * matCall(L_in,4,site_index,jump0);
}

matrix_complex simulation::georgiGlashowActionLambdaPart(long unsigned int site_index, const lattice& L_in) const
{
        matrix_complex phi;
        int jump0[4] = {0};
        phi =  matCall(L_in,4,site_index,jump0);
        return 2.0 * lambda * phi *(phi*phi).trace();
}
