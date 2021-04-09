/*  This code contains all my routines for doing Metropolis Hastings.
   My intention is to mainly use MHMC just to confirm the results of my HMC code.
   -Paul  */


#include <iostream>
#include <algorithm>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h>
#include "sim.h"
#include "stopwatch.h"
#include "lattice.h"

//=============================================================================
//The two following routines evolve either one link field or one Higgs field,
//depending on the presence of the second argument.
void simulation::evolveFieldMHMC(long unsigned int site_index, int dir)
{
        Ltemp[0].site[site_index].link[dir] = smallSU2Matrix(randomGenerator)*L.site[site_index].link[dir];
        std::cout << "in evolve step, site index is " << site_index <<"\n";
}

void simulation::evolveFieldMHMC(long unsigned int site_index)
{
        Ltemp[0].site[site_index].higgs = smallHermitianMatrix(randomGenerator) + L.site[site_index].higgs;
}


//=============================================================================
//The two following functions are meant to the take the difference in the
//Georgi-Glashow action after evolving just one site field. Creating a separate
//routine to do this difference is necessary because recalculating the action
//after every update would be very computationally intensive.
double simulation::actionDifference(long unsigned int site_index, int dir)
{
        std::cout << "in actiondifference. site_index is " << site_index  << "\n";
        double oldMixed, newMixed;
        double oldPlaq, newPlaq;
        double oldS, newS;
        double invg2;
        long unsigned int tempSiteIndex;
        int dir1;
        invg2 = 1.0d/(g*g);

        std::cout << "in action diff: post initialization" << std::endl;
        oldMixed = mixedGaugeHiggsTerm(L,site_index,dir).trace().real();
        newMixed = mixedGaugeHiggsTerm(Ltemp[0],site_index,dir).trace().real();
        std::cout << "after gaugfe terms" << std::endl;
        oldPlaq = 0.0;
        newPlaq = 0.0;

        FORALLDIRBUT(dir1,dir)
        {
                std::cout << "plaq" << std::endl;
                oldPlaq +=  plaquette(L,site_index, dir, dir1).trace().real();
                newPlaq +=  plaquette(Ltemp[0],site_index, dir, dir1).trace().real();
                std::cout << "2" << std::endl;
                tempSiteIndex = L.jumpIndex(site_index,dir1,-1);
                std::cout << "3" << std::endl;
                std::cout << "in actiondiff function pre plaquette. temp_site_index is " << tempSiteIndex  << "\n";
                oldPlaq +=  plaquette(L,tempSiteIndex, dir, dir1).trace().real();
                newPlaq +=  plaquette(Ltemp[0],tempSiteIndex, dir, dir1).trace().real();
        }

        oldS = -2.0*dsv*oldMixed - invg2*oldPlaq;
        newS = -2.0*dsv*newMixed - invg2*newPlaq;
        return newS - oldS;
}

double simulation::actionDifference(long unsigned int site_index)
{
        int jump[4]={0};
        matrix_complex phi, phiNew;
        double phi2T,phiNew2T;
        double mixedTerm, mixedTermNew;
        double oldS, newS;
        int dir;
        long unsigned int tempSiteIndex;

        phi = matCall(L,4,site_index,jump);
        phiNew = matCall(Ltemp[0],4,site_index,jump);
        phi2T = (phi*phi).trace().real();
        phiNew2T = (phiNew*phiNew).trace().real();

        mixedTerm = 0.0d; mixedTermNew = 0.0d;
        FORALLDIR(dir)
        {
                mixedTerm += mixedGaugeHiggsTerm(L,site_index,dir).trace().real();
                mixedTermNew += mixedGaugeHiggsTerm(Ltemp[0],site_index,dir).trace().real();
                tempSiteIndex = L.jumpIndex(site_index,dir,-1);
                mixedTerm += mixedGaugeHiggsTerm(L,tempSiteIndex,dir).trace().real();
                mixedTermNew += mixedGaugeHiggsTerm(Ltemp[0],tempSiteIndex,dir).trace().real();
        }
        oldS = (8+m2)* phi2T + lambda * phi2T * phi2T - 2.0*mixedTerm;
        newS = (8+m2)* phiNew2T + lambda * phiNew2T * phiNew2T - 2.0*mixedTermNew;
        return newS - oldS;
}

//=============================================================================
//The following two routines determines whether or not the proposed moved should be accepted.
void simulation::acceptOrReject(long unsigned int site_index,int dir)
{
        std::cout << "in accept function. site_index is " << site_index  << "\n";
        double deltaS = actionDifference(site_index,dir); //HERE

        double Prb = std::min( std::exp(-deltaS),1.0   );
        double Rnd = uniformReal(randomGenerator,0.0,1.0);

        if(Rnd < Prb) //Accept
        {
                AcceptanceCounter(true);
                L.site[site_index].link[dir] = Ltemp[0].site[site_index].link[dir];
        }
        else //Reject
        {
                AcceptanceCounter(false);
                Ltemp[0].site[site_index].link[dir] = L.site[site_index].link[dir];
        }
}

void simulation::acceptOrReject(long unsigned int site_index)
{
        std::cout << "in accept function. site_index is " << site_index  << "\n";
        double deltaS = actionDifference(site_index);
        double Prb = std::min( std::exp(-deltaS),1.0   );
        double Rnd = uniformReal(randomGenerator,0.0,1.0);
        if(Rnd < Prb) //Accept
        {
                AcceptanceCounter(true);
                L.site[site_index].higgs = Ltemp[0].site[site_index].higgs;
        }
        else //Reject
        {
                AcceptanceCounter(false);
                Ltemp[0].site[site_index].higgs = L.site[site_index].higgs;
        }
}
//=============================================================================
//The following routine makes one pass through every field on every site on the
//lattice. It tries to evolve each site and then moves to the next. Notice that
//Higgs/gauge field evolution can be done separately.
void simulation::sweepMHMC()
{
        unsigned long int site_index;
        int dir;
        for(site_index = 0; site_index< L.nsites; site_index++)
        {

                #ifdef __GAUGE_EVOLUTION__ //Turns gauge evolution on/off
                FORALLDIR(dir)
                {

                        evolveFieldMHMC(site_index, dir);

                        acceptOrReject(site_index,dir);


                }
                #endif

                #ifdef __HIGGS_EVOLUTION__ //Turns Higgs evolution on/off
                evolveFieldMHMC(site_index);
                acceptOrReject(site_index);
                #endif
                #ifdef __CHECK_LATTICE__ //If defined, checks to make all fields work as expected.
                isLatticeConsistent(L);
                #endif
        }
}

//=============================================================================
//This routine performs multiple MHMC sweeps through the lattice.
void simulation::multiSweepMHMC(int Nsweeps)
{
        for(int i = 0; i < Nsweeps; i++)
        {
                sweepMHMC();

        }

        //printAcceptance();
}
