#pragma once

#ifndef __SIMULATION__

#define __SIMULATION__
//#define __GAUGE_EVOLUTION__
#define __HIGGS_EVOLUTION__

#define __CHECK_NAN__
#define __CHECK_SU2__
#define __CHECK_TRACELESS__
#define __CHECK_HERMITIAN__
#define __CHECK_LATTICE__

#define CLOSETOZERO 1.0E-6

#include <random>
#include "lattice.h"
#include "rand.h"
#define DEFAULT_STEPS 100
#define DEFAULT_STEP_SIZE 1.0/DEFAULT_STEPS

#define DEFAULT_LAMBDA  0.0 //0.1
#define DEFAULT_M2 0.0  //-0.2
#define DEFAULT_STARTING_G 0.4472135955 //1.318 //
#define matCall ((*this).*boundary_condition)
typedef  const matrix_complex (simulation::*simMatrixReturn)(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4])const; //matrix_num represents either link variable number or (5) the higgs field

class simulation {
public:
//Initialization
simulation();
simulation(const simulation& sim);
simulation(const simulation& sim, char boundaryType);
simulation(double m2_in,double lambda_in,double g_in);
simulation(const lattice& L_in);
simulation(double m2_in,double lambda_in,double g_in, const lattice& L_in);
simulation(double m2_in,double lambda_in,double g_in, const lattice& L_in, char boundaryType);
~simulation();
//HMC routines
void updateFields(long unsigned site_index,double time_step, const Plattice & P_in,const lattice & L_in, lattice& L_out );
void updateMomenta(long unsigned site_index,double time_step, const Plattice & P_in,const lattice & L_in, Plattice& P_out );
void initializeHMC();
double runLeapfrogSimulation();
void copyLatticeAndRefresh(const lattice & L_in,lattice & L_out);
void copyLatticePlattice(const lattice & L_in, const Plattice & P_in, lattice & L_out, Plattice & P_out);
void leapfrogOneStep();
void wholeStepEvolve(lattice L_in, Plattice P_in, lattice L_out, Plattice P_out);
bool metropolisDecision();
//Action functions
double georgiGlashowLagrangianDensity(long unsigned int) const;
double georgiGlashowLagrangianDensity(long unsigned int, const lattice& L_in) const;
double georgiGlashowAction() const;
double georgiGlashowAction(const lattice& L_in) const;
double kineticTerm(const Plattice& P_in) const;
double georgiGlashowHamiltonian(const lattice& L_in, const Plattice& P_in) const;
double Hamiltonian() const;
const matrix_complex georgiGlashowActionLinkDerivative(long unsigned int site_index, int dir, const lattice& L_in) const;
const matrix_complex georgiGlashowActionPhiDerivative(long unsigned int site_index, const lattice& L_in) const;
//Observables
double averagePlaquettes() const;   //
const matrix_complex averagePhi() const;
const matrix_complex averagePhi2() const;
//Setup functions
void setupBoundaryConditions();
void setupBoundaryConditions( char boundaryType);
void setupParams(double m2_in,double lambda_in, double g_in);
void setupParams(double m2_in);
void setupSteps(int Nsteps);
void resetMomenta();
void resetAcceptanceCounter();
const matrix_complex periodicBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]) const;
//const matrix_complex cBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]);
//const matrix_complex twistedBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]);
void AcceptanceCounter(bool updateStatus);
void printAcceptance() const;
void printSite(long unsigned int site_index) const;
void printDerivatives(long unsigned int site_index) const;
//private:
int steps;
int nAccepts, nRejects;
double stepSize;
double m2, lambda, g;
const matrix_complex plaquette(long unsigned site_index, int dir1, int dir2) const;
const matrix_complex plaquette(const lattice& L_in,long unsigned site_index, int dir1, int dir2) const;
lattice L, Ltemp[2];
Plattice P, Ptemp[2];
std::mt19937_64 randomGenerator;
simMatrixReturn boundary_condition;
};

matrix_complex CayleyHamiltonExp(const matrix_complex & A);
bool isSU2(const matrix_complex & A);
bool isTraceless(const matrix_complex & A);
bool isHermitian(const matrix_complex &A);
bool isLatticeConsistent(const lattice &L_in );
#endif
