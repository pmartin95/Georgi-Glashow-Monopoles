#include<iostream>
#include<algorithm>
#include <unsupported/Eigen/MatrixFunctions>
#include "sim.h"

//Initialization
simulation::simulation()
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
  L = lattice(randomGenerator);
  Lcopy = L;
  startMomentum = Plattice(randomGenerator);
  endMomentum = startMomentum;
  setupBoundaryConditions('p');
}

simulation::simulation(const simulation& sim)
{
  nAccepts = 0;
  nRejects = 0;
  randomGenerator = sim.randomGenerator;
  steps = sim.steps;
  stepSize = sim.stepSize;
  m2 = sim.m2;
  lambda = sim.lambda;
  g = sim.g;
  L = sim.L;
  Lcopy = L;
  startMomentum = sim.startMomentum;
  endMomentum = sim.endMomentum; //This might need to be changed later
  boundary_condition = sim.boundary_condition;
}

simulation::simulation(const simulation& sim, char boundaryType)
{
  nAccepts = 0;
  nRejects = 0;
  randomGenerator = sim.randomGenerator;
  steps = sim.steps;
  stepSize = sim.stepSize;
  m2 = sim.m2;
  lambda = sim.lambda;
  g = sim.g;
  L = sim.L;
  Lcopy = L;
  startMomentum = sim.startMomentum;
  endMomentum = sim.endMomentum; //This might need to be changed later
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
  Lcopy = L;
  startMomentum = Plattice(randomGenerator);
  endMomentum = startMomentum;
  setupBoundaryConditions('p');
}
simulation::simulation(double m2_in,double lambda_in,double g_in, const lattice& L_in)
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
  L = L_in;
  Lcopy = L;
  startMomentum = Plattice(randomGenerator);
  endMomentum = startMomentum;
  setupBoundaryConditions('p');
}
simulation::simulation(double m2_in,double lambda_in,double g_in, const lattice& L_in, char boundaryType)
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
  L = L_in;
  Lcopy = L;
  startMomentum = Plattice(randomGenerator);
  endMomentum = startMomentum;
  setupBoundaryConditions(boundaryType);
}

//HMC routines
void simulation::runLeapfrogSimulation()
{
  for(int i;i<steps;i++)
    leapfrogOneStep();
  if(metropolisDecision())
  {
    L = Lcopy;
    resetMomenta();
  }
  else
  {
    Lcopy = L;
    resetMomenta();
  }
}

void simulation::leapfrogOneStep()
{

  long unsigned int site_index;
  int dir;
  for(site_index=0; site_index < L.nsites;site_index++)
  {
    FORALLDIR(dir)
    {
      endMomentum.site[site_index].link[dir] = endMomentum.site[site_index].link[dir] - stepSize/2.0 * georgiGlashowActionLinkDerivative(site_index, dir, Lcopy) ;
      Lcopy.site[site_index].link[dir] = Lcopy.site[site_index].link[dir] + stepSize*endMomentum.site[site_index].link[dir] ;
    }
    endMomentum.site[site_index].higgs = endMomentum.site[site_index].higgs - stepSize/2.0 * georgiGlashowActionPhiDerivative(site_index, Lcopy);
    Lcopy.site[site_index].higgs = Lcopy.site[site_index].higgs + stepSize * endMomentum.site[site_index].higgs;
  }
  for(site_index=0; site_index < L.nsites;site_index++)
  {
    FORALLDIR(dir)
    {
      endMomentum.site[site_index].link[dir] = endMomentum.site[site_index].link[dir] - stepSize/2.0 * georgiGlashowActionLinkDerivative(site_index, dir, Lcopy);
    }
    endMomentum.site[site_index].higgs = endMomentum.site[site_index].higgs - stepSize/2.0 * georgiGlashowActionPhiDerivative(site_index, Lcopy);
  }
}

bool simulation::metropolisDecision()
{
  double expResult;
  double randomDecider;
  //should return true if to accept configuration, otherwise replace with copy and start over
  expResult = std::min(exp((georgiGlashowHamiltonian(L,endMomentum) - georgiGlashowHamiltonian(Lcopy,startMomentum))), 1.0);
  randomDecider = uniformReal(randomGenerator);
  if(expResult > randomDecider)
  {
    nAccepts++;
    return true;
  }
  else
  {
    nRejects++;
    return false;
  }
}

//Action functions
double simulation::georgiGlashowLagrangianDensity(long unsigned int site_index) const
{
  double total = 0.0, subtotal = 0.0;
  int i,j;
  int jumpNone[4] = {0,0,0,0};
  long unsigned int temp_index;
  matrix_complex phi  = ((*this).*boundary_condition)(L,4,site_index,jumpNone);
  matrix_complex phi_temp;
  matrix_complex link1,link2, link3, link4;

  total += 8.0 * (phi * phi).trace().real();
  FORALLDIR(i)
  {
    int jump1[4] = {0};
    jump1[i]+=1;
    phi_temp = ((*this).*boundary_condition)(L,4, site_index, jump1);
    link1 = ((*this).*boundary_condition)(L,i,site_index,jumpNone);
    total -= 2.0 * (phi * link1 * phi_temp * link1.adjoint() ).trace().real();
  }
  FORALLDIRLESSTHAN(i,j)
  {
    subtotal += 2.0 - (plaquette(site_index, i,j) ).trace().real();
  }
  total += subtotal * 2.0 /(g*g);
  total += m2 *(phi * phi).trace().real();
  subtotal = (phi * phi).trace().real();
  total += lambda * subtotal * subtotal;
  return total;
}

double simulation::georgiGlashowLagrangianDensity(long unsigned int site_index, const lattice& L_in) const
{
  double total = 0.0, subtotal = 0.0;
  int i,j;
  int jumpNone[4] = {0,0,0,0};
  long unsigned int temp_index;
  matrix_complex phi  = ((*this).*boundary_condition)(L_in,4,site_index,jumpNone);
  matrix_complex phi_temp;
  matrix_complex link1,link2, link3, link4;

  total += 8.0 * (phi * phi).trace().real();
  FORALLDIR(i)
  {
    int jump1[4] = {0};
    jump1[i]+=1;
    phi_temp = ((*this).*boundary_condition)(L_in,4, site_index, jump1);
    link1 = ((*this).*boundary_condition)(L_in,i,site_index,jumpNone);
    total -= 2.0 * (phi * link1 * phi_temp * link1.adjoint() ).trace().real();
  }
  FORALLDIRLESSTHAN(i,j)
  {
    subtotal += 2.0 - (plaquette(site_index, i,j) ).trace().real();
  }
  total += subtotal * 2.0 /(g*g);
  total += m2 *(phi * phi).trace().real();
  subtotal = (phi * phi).trace().real();
  total += lambda * subtotal * subtotal;
  return total;
}

double simulation::georgiGlashowAction() const
{
  double total = 0.0;
  for(long unsigned int i = 0; i<L.nsites;i++)
    total += georgiGlashowLagrangianDensity(i);
  return total;
}
double simulation::georgiGlashowAction(const lattice& L_in) const
{
  double total = 0.0;
  for(long unsigned int i = 0; i<L_in.nsites;i++)
    total += georgiGlashowLagrangianDensity(i,L_in);
  return total;
}

double simulation::georgiGlashowHamiltonian(lattice& L_in, Plattice& P_in) const
{
  matrix_complex momenta_matrix;
  double field_total, momenta_total;
  long unsigned int site_index;
  int i;

  momenta_total = 0.0;
  momenta_matrix.setZero();

  field_total = georgiGlashowAction(L_in);

  for(site_index=0;site_index < L_in.nsites; site_index++)
  {
    FORALLDIR(i)
      momenta_matrix += P_in.site[site_index].link[i].transpose() *  P_in.site[site_index].link[i];
    momenta_matrix +=   P_in.site[site_index].higgs.transpose() *  P_in.site[site_index].higgs;
  }
  momenta_total = momenta_matrix.trace().real()/2.0;
  return field_total + momenta_total;
}

const matrix_complex simulation::georgiGlashowActionLinkDerivative(long unsigned int site_index, int dir, const lattice& L_in) const //
{
  matrix_complex temp1,temp2,temp3,temp4,Iden;
  matrix_complex subtotal1,subtotal2;
  int jump1[4] = {0}, jump2[4] = {0}, jump3[4] = {0}, jump4[4] = {0};
  int jumpNone[4] = {0,0,0,0};
  int temp_dir;
  complex<double> traced_part;
  Iden.setIdentity();
  jump1[dir]++;
  FORALLDIRBUT(temp_dir,dir)
  {
      jump2[temp_dir]++;
      jump3[temp_dir]--;
      jump4[dir]++;
      jump4[temp_dir]--;
      temp1 = ((*this).*boundary_condition)(L_in,dir,site_index,jumpNone) * ((*this).*boundary_condition)(L_in,temp_dir,site_index,jump1) * ((*this).*boundary_condition)(L_in,dir,site_index,jump2).adjoint() * ((*this).*boundary_condition)(L_in,temp_dir,site_index,jumpNone).adjoint();
      temp2 = ((*this).*boundary_condition)(L_in,temp_dir,site_index,jumpNone) * ((*this).*boundary_condition)(L_in,dir,site_index,jump2) * ((*this).*boundary_condition)(L_in,temp_dir,site_index,jump1).adjoint() * ((*this).*boundary_condition)(L_in,dir,site_index,jumpNone).adjoint();
      temp3 = ((*this).*boundary_condition)(L_in,temp_dir,site_index,jump3).adjoint() * ((*this).*boundary_condition)(L_in,dir,site_index,jump3) * ((*this).*boundary_condition)(L_in,temp_dir,site_index,jump4) * ((*this).*boundary_condition)(L_in,dir,site_index,jumpNone).adjoint();
      temp4 = ((*this).*boundary_condition)(L_in,dir,site_index,jumpNone) * ((*this).*boundary_condition)(L_in,temp_dir,site_index,jump4).adjoint() * ((*this).*boundary_condition)(L_in,dir,site_index,jump3).adjoint() * ((*this).*boundary_condition)(L_in,temp_dir,site_index,jump3);
      std::fill(std::begin(jump2),std::end(jump2),0);
      std::fill(std::begin(jump3),std::end(jump3),0);
      std::fill(std::begin(jump4),std::end(jump4),0);
  }
  traced_part = (temp1 - temp2 - temp3 + temp4).trace();
  subtotal1 = 2.0d*(temp1 - temp2 - temp3 + temp4) - traced_part * Iden;

  temp1 = ((*this).*boundary_condition)(L_in,dir,site_index,jumpNone) * ((*this).*boundary_condition)(L_in,4,site_index,jump1) * ((*this).*boundary_condition)(L_in,dir,site_index,jumpNone).adjoint()* ((*this).*boundary_condition)(L_in,4,site_index,jumpNone);
  temp2 = ((*this).*boundary_condition)(L_in,4,site_index,jumpNone) * ((*this).*boundary_condition)(L_in,dir,site_index,jumpNone) * ((*this).*boundary_condition)(L_in,4,site_index,jump1)* ((*this).*boundary_condition)(L_in,dir,site_index,jumpNone).adjoint();
  subtotal2 = temp1 - temp2;
  return subtotal1 * complex<double>(0.0,1.0/(g*g)) + subtotal2 * complex<double>(0.0,0.5);
}


const matrix_complex simulation::georgiGlashowActionPhiDerivative(long unsigned int site_index, const lattice& L_in) const
{
    matrix_complex temp1, temp2, temp3,Iden;
    int temp_dir;
    int jumpNone[4] = {0,0,0,0};
    double traced_part;
    temp1.setZero();
    temp2.setZero();
    temp3.setZero();
    Iden.setIdentity();

    FORALLDIR(temp_dir)
    {
      int jump1[4] = {0}, jump2[4] ={0};
      jump1[temp_dir]++;
      jump2[temp_dir]--;
      temp2 += ((*this).*boundary_condition)(L_in,temp_dir,site_index,jumpNone) *  ((*this).*boundary_condition)(L_in,4,site_index,jump1) * ((*this).*boundary_condition)(L_in,temp_dir,site_index,jumpNone).adjoint();
      temp3 += ((*this).*boundary_condition)(L_in,temp_dir,site_index,jump2).adjoint() *  ((*this).*boundary_condition)(L_in,4,site_index,jump2)* ((*this).*boundary_condition)(L_in,temp_dir,site_index,jump2);
    }
    temp1 +=8.0*lambda* ((*this).*boundary_condition)(L_in,4,site_index,jumpNone);
    temp1 = (temp1*(temp1*temp1).trace() + 4*m2*temp1+ 32*temp1  );
    traced_part = 2.0 * (temp2 + temp3).trace().real();
    return temp1 - 4*(temp2 + temp3) + traced_part * Iden;
}

//Observables
double simulation::averagePlaquettes() const
{
  int dir1,dir2;
  double subtotal = 0.0d;
  unsigned long int site_index;
  FORALLDIRLESSTHAN(dir1,dir2)
    subtotal += plaquette(site_index,dir1,dir2).trace().real();
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

void simulation::resetMomenta()
{
  startMomentum = Plattice(randomGenerator);
  endMomentum = startMomentum;
}


const matrix_complex simulation::periodicBoundaryCondition(const lattice& L_in,int matrix_num, unsigned long int index,const int jump[4]) const
{
  int i;
  bool shouldBreak;
  FORALLDIR(i)
  {
    shouldBreak = jump[i] != 0;
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
  if(x[dir] >= L_in.ns[dir])
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
void simulation::printAcceptance() const
{
  if(nAccepts + nRejects > 0)
    std::cout << static_cast<float>(nAccepts)/ static_cast<float>(nAccepts+nRejects) <<std::endl;
  else
    std::cout << "There have been no acceptances or rejections.\n";
}

const matrix_complex simulation::plaquette(long unsigned site_index, int dir1, int dir2) const
{
  matrix_complex u1, u2, u3, u4;
  int jump1[4] ={0}, jump2[4]={0}, jump3[4]={0}, jump4[4]={0};
  jump1[dir1] = 0;
  u1 = ((*this).*boundary_condition)(L,dir1, site_index,jump1);
  jump2[dir1] = 1;
  u2 = ((*this).*boundary_condition)(L,dir2, site_index,jump2);
  jump3[dir2] = 1;
  u3 = ((*this).*boundary_condition)(L,dir1, site_index,jump3).adjoint();
  jump4[dir2] = 0;
  u4 = ((*this).*boundary_condition)(L,dir2, site_index,jump4).adjoint();

  return u1*u2*u3*u4;
}

const matrix_complex simulation::plaquette(const lattice& L_in,long unsigned site_index, int dir1, int dir2) const
{
  matrix_complex u1, u2, u3, u4;
  int jump1[4] ={0}, jump2[4]={0}, jump3[4]={0}, jump4[4]={0};

  jump1[dir1] = 0;
  u1 = ((*this).*boundary_condition)(L_in,dir1, site_index,jump1);
  jump2[dir1] = 1;
  u2 = ((*this).*boundary_condition)(L_in,dir2, site_index,jump2);
  jump3[dir2] = 1;
  u3 = ((*this).*boundary_condition)(L_in,dir1, site_index,jump3).adjoint();
  jump4[dir2] = 0;
  u4 = ((*this).*boundary_condition)(L_in,dir2, site_index,jump4).adjoint();

  return u1*u2*u3*u4;
}
