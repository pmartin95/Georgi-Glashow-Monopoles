#pragma once

#ifndef __SIMULATION__

#define __SIMULATION__

#include <random>
#include "lattice.h"
#include "rand.h"
#define DEFAULT_STEPS 200
#define DEFAULT_STEP_SIZE 0.05
#define DEFAULT_LAMBDA 1.0
#define DEFAULT_M2 1.0
#define DEFAULT_STARTING_G 1.0
#define matCall ((*this).*boundary_condition) 
typedef  const matrix_complex (simulation::*simMatrixReturn)(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4])const; //matrix_num represents either link variable number or (5) the higgs field

class simulation{
public:
  //Initialization
  simulation();
  simulation(const simulation& sim);
  simulation(const simulation& sim, char boundaryType);
  simulation(double m2_in,double lambda_in,double g_in);
  simulation(double m2_in,double lambda_in,double g_in, const lattice& L_in);
  simulation(double m2_in,double lambda_in,double g_in, const lattice& L_in, char boundaryType);
  //HMC routines
  void runLeapfrogSimulation();
  void leapfrogOneStep();
  bool metropolisDecision();
  //Action functions
  double georgiGlashowLagrangianDensity(long unsigned int) const;
  double georgiGlashowLagrangianDensity(long unsigned int, const lattice& L_in) const;
  double georgiGlashowAction() const;
  double georgiGlashowAction(const lattice& L_in) const;
  double georgiGlashowHamiltonian(lattice& L_in, Plattice& P_in) const;
  const matrix_complex georgiGlashowActionLinkDerivative(long unsigned int site_index, int dir, const lattice& L_in) const;
  const matrix_complex georgiGlashowActionPhiDerivative(long unsigned int index, const lattice& L_in) const;
  //Observables
  double averagePlaquettes() const;
  //Setup functions
  void setupBoundaryConditions();
  void setupBoundaryConditions( char boundaryType);
  void setupParams(double m2_in,double lambda_in, double g_in);
  void setupParams(double m2_in);
  void resetMomenta();//
  const matrix_complex periodicBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]) const;
  //const matrix_complex cBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]);
  //const matrix_complex twistedBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]);
  void printAcceptance() const;
private:
  int steps;
  int nAccepts, nRejects;
  double stepSize;
  double m2, lambda, g;
  const matrix_complex plaquette(long unsigned index, int dir1, int dir2) const;
  const matrix_complex plaquette(const lattice& L_in,long unsigned index, int dir1, int dir2) const;
  lattice L, Lcopy;
  Plattice startMomentum, endMomentum;
  std::mt19937_64 randomGenerator;
  simMatrixReturn boundary_condition;
};

#endif
