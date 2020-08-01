#include "lattice.h"
#include "rand.h"

#ifndef __SIMULATION__
#define __SIMULATION__

#define DEFAULT_STEPS 200
#define DEFAULT_STEP_SIZE 0.05
#define DEFAULT_LAMBDA 1.0
#define DEFAULT_M2 1.0
#define DEFAULT_STARTING_G 1.0

typedef  const matrix_complex (simulation::*simMatrixReturn)(int matrix_num, unsigned long int index, int dir, int jump); //matrix_num represents either link variable number or (5) the higgs field

class simulation{
public:
  //Initialization
  simulation(); //
  simulation(const simulation& sim); //
  simulation(const simulation& sim, char boundaryType); //
  simulation(double m2_in,double lambda_in,double g_in); //
  simulation(double m2_in,double lambda_in,double g_in, const lattice& L_in); //
  simulation(double m2_in,double lambda_in,double g_in, const lattice& L_in, char boundaryType); //
  //HMC routines
  void runLeapfrogSimulation(); //
  void leapfrogOneStep(int i); //
  void evolveMomentum(int i);//
  void evolveFields(int i);//
  bool metropolisDecision();//
  //Action functions
  double georgiGlashowLagrangianDensity(long unsigned int) const;
  double georgiGlashowAction() const;
  double georgiGlashowHamiltonian(lattice * L_in, Plattice *P_in) const;
  const matrix_complex georgiGlashowActionLinkDerivative(long unsigned int) const; //
  const matrix_complex georgiGlashowActionPhiDerivative(long unsigned int) const; //
  //Observables
  double averagePlaquettes() const;//
  //Setup functions
  void setupBoundaryConditions();
  void setupBoundaryConditions( char boundaryType);
  void setupParams(int m2,int lambda, int g); //
  const matrix_complex periodicBoundaryCondition(int matrix_num, unsigned long int index, int dir, int jump);
  //const matrix_complex cBoundaryCondition(int matrix_num, unsigned long int index, int dir, int jump);
  //const matrix_complex twistedBoundaryCondition(int matrix_num, unsigned long int index, int dir, int jump);
private:
  int steps;
  int nAccepts, nRejects;
  double stepSize;
  double m2, lambda, g;
  const matrix_complex plaquette(long unsigned index, int dir1, int dir2) const;
  lattice L, Lcopy;
  Plattice startMomentum, endMomentum;
  URNG randomGenerator;
  simMatrixReturn boundary_condition;
};

#endif
