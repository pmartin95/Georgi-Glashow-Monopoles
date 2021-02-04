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


void simulation::evolveFieldMHMC(long unsigned int site_index, int dir)
{
        Ltemp[0].site[site_index].link[dir] = smallSU2Matrix(randomGenerator)*L.site[site_index].link[dir];
}
void simulation::evolveFieldMHMC(long unsigned int site_index)
{
        Ltemp[0].site[site_index].higgs = smallHermitianMatrix(randomGenerator) + L.site[site_index].higgs;
}
double simulation::actionDifference(long unsigned int site_index, int dir)
{
        double oldMixed, newMixed;
        double oldPlaq, newPlaq;
        double oldS, newS;
        double invg2;
        long unsigned int tempSiteIndex;
        int dir1;
        invg2 = 1.0d/(g*g);
        oldMixed = mixedGaugeHiggsTerm(L,site_index,dir).trace().real();
        newMixed = mixedGaugeHiggsTerm(Ltemp[0],site_index,dir).trace().real();

        oldPlaq = 0.0;
        newPlaq = 0.0;
        FORALLDIRBUT(dir1,dir)
        {
                oldPlaq +=  plaquette(L,site_index, dir, dir1).trace().real();
                newPlaq +=  plaquette(Ltemp[0],site_index, dir, dir1).trace().real();
                tempSiteIndex = L.jumpIndex(site_index,dir1,-1);
                oldPlaq +=  plaquette(L,tempSiteIndex, dir, dir1).trace().real();
                newPlaq +=  plaquette(Ltemp[0],tempSiteIndex, dir, dir1).trace().real();
        }

        oldS = -2.0*oldMixed - invg2*oldPlaq;
        newS = -2.0*newMixed - invg2*newPlaq;
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

void simulation::acceptOrReject(long unsigned int site_index,int dir)
{
        double deltaS = actionDifference(site_index,dir);
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

void simulation::sweepMHMC()
{
        unsigned long int site_index;
        int dir;
        for(site_index = 0; site_index< L.nsites; site_index++)
        {
                FORALLDIR(dir)
                {
                        evolveFieldMHMC(site_index, dir);
                        acceptOrReject(site_index,dir);
                }
                evolveFieldMHMC(site_index, dir);
                acceptOrReject(site_index,dir);
        }
}

void simulation::multiSweepMHMC(int Nsweeps)
{
        for(int i = 0; i < Nsweeps; i++)
                sweepMHMC();
        printAcceptance();
}
