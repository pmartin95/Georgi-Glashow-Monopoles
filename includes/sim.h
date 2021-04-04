#pragma once

#ifndef __SIMULATION__

#define __SIMULATION__
#define __GAUGE_EVOLUTION__
//#define __HIGGS_EVOLUTION__
//#define __HIGGS_GAUGE_MIXED_TERM__

// #define __CHECK_NAN__
// #define __CHECK_SU2__
// #define __CHECK_TRACELESS__
// #define __CHECK_HERMITIAN__
// #define __CHECK_LATTICE__

#define CLOSETOZERO 1.0E-6

#include <random>
#include <string>
#include "lattice.h"
#include "rand.h"
#include <iomanip>
#define DEFAULT_STEPS 10
#define DEFAULT_STEP_SIZE 1.0/DEFAULT_STEPS

#define FIRST_TERM_PARAM 1.0 //for turning on an off the first part of the action
#define DEFAULT_LAMBDA  0.1
#define DEFAULT_M2 -0.2
#define DEFAULT_STARTING_G 0.4472135955 //1.318 //


#define matCall ((*this).*boundary_condition)

typedef double (*oneVarDoubleFunc) (double);

typedef  const matrix_complex (simulation::*simMatrixReturn)(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4])const; //matrix_num represents either link variable number or (5) the higgs field

typedef struct data_point
{
        double g_value;
        double m2_value;
        double lambda_value;
        double average_plaquette_value;
        double average_phi2_value;
} data_point_t;

typedef struct schedule_element
{
        double g_value;
        double m2_value;
        double lambda_value;
} schedule_element_t;

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
//Metropolis Hastings Routines
void evolveFieldMHMC(long unsigned int site_index, int dir);
void evolveFieldMHMC(long unsigned int site_index);
double actionDifference(long unsigned int site_index, int dir);
double actionDifference(long unsigned int site_index);
void acceptOrReject(long unsigned int site_index,int dir);
void acceptOrReject(long unsigned int site_index);
void sweepMHMC();
void multiSweepMHMC(int Nsweeps);
//Action functions
double georgiGlashowLagrangianDensity(long unsigned int) const;
double georgiGlashowLagrangianDensity(long unsigned int, const lattice& L_in) const;
double georgiGlashowAction() const;
double georgiGlashowAction(const lattice& L_in) const;
double kineticTerm(const Plattice& P_in) const;
double georgiGlashowHamiltonian(const lattice& L_in, const Plattice& P_in) const;
double Hamiltonian() const;
const matrix_complex georgiGlashowActionLinkDerivative(long unsigned int site_index, int dir, const lattice& L_in) const;
const matrix_complex georgiGlashowActionPureGaugeDerivative(long unsigned int site_index, int dir, const lattice& L_in) const;
const matrix_complex georgiGlashowActionMixedGaugeDerivative(long unsigned int site_index, int dir, const lattice& L_in) const;
const matrix_complex georgiGlashowActionPhiDerivative(long unsigned int site_index, const lattice& L_in) const;
const matrix_complex georgiGlashowActionPhiKineticPart(long unsigned int site_index, const lattice& L_in) const;
const matrix_complex georgiGlashowActionPhiMPart(long unsigned int site_index, const lattice& L_in) const;
const matrix_complex georgiGlashowActionLambdaPart(long unsigned int site_index, const lattice& L_in) const;
//Observables
double averagePlaquettes() const;   //
const matrix_complex averagePhi() const;
const matrix_complex averagePhi2() const;
double CreutzRatio() const;
//Setup and Reset functions
void setupBoundaryConditions();
void setupBoundaryConditions( char boundaryType);
void setupParams(double m2_in,double lambda_in, double g_in);
void setupParams(double m2_in);
void switchDSV();
void setupSteps(int Nsteps);
void resetMomenta();
void resetAcceptanceCounter();
//Boundary Conditions
const matrix_complex directMatCall(const lattice& L_in,unsigned long int site_index,int matrix_num) const;
int shiftToLattice(const lattice& L_in,int coordinate, int dir) const;
int incrementCoordinate(const lattice& L_in,int coordinate,int dir) const;
const matrix_complex periodicBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]) const;
const matrix_complex cBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]) const;
const matrix_complex cTwist(const matrix_complex& mat_A, int matrix_num) const;
const matrix_complex cTwistGaugeField(const matrix_complex& mat_A) const;
const matrix_complex cTwistHiggsField(const matrix_complex& mat_A) const;
const matrix_complex twistedBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index, const int jump[4]) const;
const matrix_complex TwistField(const matrix_complex& mat_A,int matrix_number,int dir) const;
const matrix_complex xTwist(const matrix_complex& mat_A,int matrix_number) const;
const matrix_complex yTwist(const matrix_complex& mat_A,int matrix_number) const;
const matrix_complex zTwist(const matrix_complex& mat_A,int matrix_number) const;
void AcceptanceCounter(bool updateStatus);
void printAcceptance() const;
void printSite(long unsigned int site_index) const;
void printDerivatives(long unsigned int site_index) const;
//Data Collection and Scheduling Routines
void appendDataPoint();
void resetDataPoints();
void printDataFile(const std::string& filename) const;
void inputScheduleParameters(const std::string& filename);
void reverseSchedule();
void generateScheduleFile(const std::string& filename,const std::vector<double>& gs,const std::vector<double>& m2s,const std::vector<double>& lambdas);
void loadScheduleValues(int i);
void runHMCSimulationSchedule(int init_iter,int iter, int iter_measure);
void runMHMCSimulationSchedule(int iter, int iter_measure);
//private:
int steps;
int nAccepts, nRejects;
double stepSize;
double m2, lambda, g,invg,dsv;
const matrix_complex plaquette(long unsigned site_index, int dir1, int dir2) const;
const matrix_complex plaquette(const lattice& L_in,long unsigned site_index, int dir1, int dir2) const;
const matrix_complex mixedGaugeHiggsTerm(const lattice& L_in, long unsigned site_index, int dir) const;
lattice L, Ltemp[2];
Plattice P, Ptemp[2];
std::mt19937_64 randomGenerator;
simMatrixReturn boundary_condition;
std::vector<data_point_t> data;
std::vector<schedule_element_t> schedule;

};


matrix_complex CayleyHamiltonExp(const matrix_complex & A);
bool isSU2(const matrix_complex & A);
bool isTraceless(const matrix_complex & A);
bool isHermitian(const matrix_complex &A);
bool isLatticeConsistent(const lattice &L_in );
vector<double> createDistributionVector(double a, double b, int N);
vector<double> createDistributionVector(double a, double b, int N, oneVarDoubleFunc func);
double invSquareRoot(double x);
double averageDoubleVector(std::vector<double> &V );
#endif
