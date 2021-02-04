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
        nAccepts = 0;
        nRejects = 0;
        std::mt19937_64 randTemp(seedGen());
        randomGenerator = randTemp;
        steps = DEFAULT_STEPS;
        stepSize = DEFAULT_STEP_SIZE;
        m2 = m2_in;
        lambda = lambda_in;
        g = g_in;
        L = lattice(randomGenerator);
        Ltemp[0] = L;
        Ltemp[1] = Ltemp[0];
        P = Plattice(randomGenerator);
        Ptemp[0] = P;
        Ptemp[1] = Ptemp[0];

        setupBoundaryConditions('p');
}

simulation::simulation(const lattice& L_in)
{
        nAccepts = 0;
        nRejects = 0;
        std::mt19937_64 randTemp(seedGen());
        randomGenerator = randTemp;
        steps = DEFAULT_STEPS;
        stepSize = DEFAULT_STEP_SIZE;
        m2 = DEFAULT_M2;
        lambda = DEFAULT_LAMBDA;
        g = DEFAULT_STARTING_G;
        L = L_in;
        Ltemp[0] = L;
        Ltemp[1] = Ltemp[0];
        P = Plattice(randomGenerator);
        Ptemp[0] = P;
        Ptemp[1] = Ptemp[0];
        setupBoundaryConditions('p');
}

simulation::simulation(double m2_in,double lambda_in,double g_in, const lattice& L_in)
{
        (*this) = simulation(L_in);
        m2 = m2_in;
        lambda = lambda_in;
        g = g_in;
        setupBoundaryConditions('p');
}
simulation::simulation(double m2_in,double lambda_in,double g_in, const lattice& L_in, char boundaryType)
{
        (*this) = simulation( m2_in,lambda_in,g_in, L_in  );
        setupBoundaryConditions(boundaryType);
}

simulation::~simulation()
{

}


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
        phiDerivativePart *= 2.0;
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
        phiDerivativePart = 0.0; /////

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


//Observables
double simulation::averagePlaquettes() const
{
        int dir1,dir2;
        double subtotal = 0.0;
        unsigned long int site_index;
  #pragma omp parallel private(dir1,dir2) default(shared)
        {
                double local_subtotal = 0.0;
    #pragma omp for
                for(long unsigned int i = 0; i< L.nsites; i++)
                {
                        FORALLDIRLESSTHAN(dir1,dir2)
                        local_subtotal+= plaquette(site_index,dir1,dir2).trace().real();
                }
    #pragma omp critical
                subtotal+=local_subtotal;
        }
        return subtotal / static_cast<double>(L.nsites);
}

const matrix_complex simulation::averagePhi() const
{
        matrix_complex subtotal;
        subtotal.setZero();
        unsigned long int site_index;
  #pragma omp parallel default(shared)
        {
                matrix_complex local_subtotal;
                subtotal.setZero();
    #pragma omp for
                for(long unsigned int i = 0; i< L.nsites; i++)
                {
                        local_subtotal+= L.site[i].higgs;
                }
    #pragma omp critical
                subtotal+=local_subtotal;
        }
        return subtotal / static_cast<double>(L.nsites);
}

const matrix_complex simulation::averagePhi2() const
{
        matrix_complex subtotal;
        subtotal.setZero();
        unsigned long int site_index;
  #pragma omp parallel default(shared)
        {
                matrix_complex local_subtotal;
                subtotal.setZero();
    #pragma omp for
                for(long unsigned int i = 0; i< L.nsites; i++)
                {
                        local_subtotal += L.site[i].higgs * L.site[i].higgs;
                }
    #pragma omp critical
                subtotal += local_subtotal;
        }
        return subtotal / static_cast<double>(L.nsites);
}

//Setup functions
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
void simulation::setupParams(double m2_in,double lambda_in, double g_in)
{
        m2 = m2_in;
        lambda = lambda_in;
        g = g_in;
}

void simulation::setupParams(double m2_in)
{
        m2 = m2_in;
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
void simulation::AcceptanceCounter(bool updateStatus)
{
        if(updateStatus)
                nAccepts++;
        else
                nRejects++;
}


const matrix_complex simulation::plaquette(long unsigned site_index, int dir1, int dir2) const
{
        matrix_complex u1, u2, u3, u4;
        int jump1[4] ={0}, jump2[4]={0}, jump3[4]={0}, jump4[4]={0};
        jump2[dir1] += 1; jump3[dir2] += 1;
        u1 = matCall(L,dir1, site_index,jump1);
        u2 = matCall(L,dir2, site_index,jump2);
        u3 = matCall(L,dir1, site_index,jump3).adjoint();
        u4 = matCall(L,dir2, site_index,jump4).adjoint();

        return u1*u2*u3*u4;
}

const matrix_complex simulation::plaquette(const lattice& L_in,long unsigned site_index, int dir1, int dir2) const
{
        matrix_complex u1, u2, u3, u4;
        int jump1[4] ={0}, jump2[4]={0}, jump3[4]={0}, jump4[4]={0};
        jump2[dir1] += 1; jump3[dir2] += 1;
        u1 = matCall(L_in,dir1, site_index,jump1);
        u2 = matCall(L_in,dir2, site_index,jump2);
        u3 = matCall(L_in,dir1, site_index,jump3).adjoint();
        u4 = matCall(L_in,dir2, site_index,jump4).adjoint();

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
