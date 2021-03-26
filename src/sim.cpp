#include <iostream>
#include <algorithm>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h>
#include "sim.h"
#include "stopwatch.h"

//Initialization
simulation::simulation()
{
        nAccepts = 0;
        nRejects = 0;
        steps = DEFAULT_STEPS;
        stepSize = DEFAULT_STEP_SIZE;
        m2 = DEFAULT_M2;
        lambda = DEFAULT_LAMBDA;
        g = DEFAULT_STARTING_G;
        invg = 1.0/g;
        dsv = FIRST_TERM_PARAM;
        std::mt19937_64 randTemp(seedGen());
        randomGenerator = randTemp;
        L = lattice(randomGenerator);
        Ltemp[0] = L;
        Ltemp[1] = Ltemp[0];
        P = Plattice(randomGenerator);
        Ptemp[0] = P;
        Ptemp[1] = Ptemp[0];
        setupBoundaryConditions('p');
}

simulation::simulation(const simulation& sim)
{
        nAccepts = sim.nAccepts;
        nRejects = sim.nRejects;
        steps = sim.steps;
        stepSize = sim.stepSize;
        m2 = sim.m2;
        lambda = sim.lambda;
        g = sim.g;
        invg = 1.0/g;

        randomGenerator = sim.randomGenerator;

        L = sim.L;
        Ltemp[0] = sim.Ltemp[0];
        Ltemp[1] = sim.Ltemp[1];
        P = sim.P;
        Ptemp[0] = sim.Ptemp[0];
        Ptemp[1] = sim.Ptemp[1];
        boundary_condition = sim.boundary_condition;
}

simulation::simulation(const simulation& sim, char boundaryType)
{
        (*this) = simulation(sim);
        setupBoundaryConditions(boundaryType);
}

simulation::simulation(double m2_in,double lambda_in,double g_in)
{
        (*this) = simulation();
        setupParams(m2_in,lambda_in,g_in);
}

simulation::simulation(const lattice& L_in)
{
        (*this) = simulation();
        L = L_in;
        Ltemp[0] = L;
        Ltemp[1] = Ltemp[0];
        P = Plattice(randomGenerator,L_in);
        Ptemp[0] = P;
        Ptemp[1] = Ptemp[0];
}

simulation::simulation(double m2_in,double lambda_in,double g_in, const lattice& L_in)
{
        (*this) = simulation(L_in);
        setupParams(m2_in,lambda_in,g_in);
}
simulation::simulation(double m2_in,double lambda_in,double g_in, const lattice& L_in, char boundaryType)
{
        (*this) = simulation( m2_in,lambda_in,g_in, L_in  );
        setupBoundaryConditions(boundaryType);
}

simulation::~simulation()
{

}



//Setup functions
void simulation::setupParams(double m2_in,double lambda_in, double g_in)
{
        m2 = m2_in;
        lambda = lambda_in;
        g = g_in;
        invg = 1.0/g;
}

void simulation::setupParams(double m2_in)
{
        m2 = m2_in;
}

void simulation::switchDSV()
{
        if(dsv < 1.0)
                dsv = 1.0;
        else
                dsv = 0.0;
}


void simulation::setupSteps(int Nsteps)
{
        steps = Nsteps;
        stepSize = 1.0/static_cast<double>(steps);
}
void simulation::resetMomenta()
{
        P = Plattice(randomGenerator);
        Ptemp[0] = P;
        Ptemp[1] = Ptemp[0];
}

void simulation::resetAcceptanceCounter()
{
        nAccepts = 0;
        nRejects = 0;
}


void simulation::AcceptanceCounter(bool updateStatus)
{
        if(updateStatus)
                nAccepts++;
        else
                nRejects++;
}
