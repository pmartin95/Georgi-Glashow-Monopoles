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
  L = L_in
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
  L = L_in
  Lcopy = L;
  startMomentum = Plattice(randomGenerator);
  endMomentum = startMomentum;
  setupBoundaryConditions(boundaryType);
}

//HMC routines
void simulation::runLeapfrogSimulation()
{
  for(int i;i<steps;i++)
    leapfrogOneStep(i);
}

void simulation::leapfrogOneStep(int i)
{
  evolveFields(i);
  evolveMomentum(i);
}

void evolveMomentum(int i);//
void evolveFields(int i);//
bool simulation::metropolisDecision()
{
  double beta = 4.0 /(g*g) //? Double check this
  double expResult;
  double randomDecider;
  //should return true if to accept configuration, otherwise replace with copy and start over
  expResult = std::min(exp(-beta * (georgiGlashowHamiltonian(L,endMomentum) - georgiGlashowHamiltonian(Lcopy,startMomentum))), 1.0)
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
double simulation::georgiGlashowLagrangianDensity(long unsigned int index) const
{
  double total = 0.0, subtotal = 0.0;
  int i,j;
  long unsigned int temp_index;
  matrix_complex phi = L.site[index].higgs;
  matrix_complex phi_temp;
  matrix_complex link1,link2, link3, link4;

  total += 8.0 * (phi * phi).trace().real();
  FORALLDIR(i)
  {
    phi_temp = boundary_condition(4, index, i, 1);
    link1 = L.site[index].link[i];
    total -= 2.0 * (phi * link1 * phi_temp * link1.adjoint() ).trace().real();
  }
  FORALLDIRLESSTHAN(i,j)
  {
    subtotal += 2.0 - (plaquette(index, i,j) ).trace().real();
  }
  total += subtotal * 2.0 /(g*g);
  total += m2 *(phi * phi).trace().real()
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

double simulation::georgiGlashowHamiltonian(lattice * L_in, Plattice *P_in) const;
const matrix_complex georgiGlashowActionLinkDerivative(long unsigned int) const; //
const matrix_complex georgiGlashowActionPhiDerivative(long unsigned int) const; //
//Observables
double simulation::averagePlaquettes() const;//
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
void simulation::setupParams(int m2_in,int lambda_in, int g_in)
{
  m2 = m2_in;
  lambda = lambda_in;
  g = g_in;
}

const matrix_complex& simulation::periodicBoundaryCondition(int matrix_num, unsigned long int index, int dir, int jump)
{
  if(jump == 0)
    return L.site[index].link[matrix_num];
  int x[4];
  long unsigned int new_index;
  L.indexToCoordinate(index, x);
  x[dir] += jump;
  if(x[dir] >= ns[dir])
    x[dir] = (x[dir] + ns[dir])%ns[dir];
  new_index = L.coordinateToIndex(x);
  if(matrix_num<4)
    return L.site[new_index].link[matrix_num];
  else
    return L.site[new_index].higgs;
}
//const matrix_complex simulation::cBoundaryCondition(int matrix_num, unsigned long int index, int dir, int jump);
//const matrix_complex simulation::twistedBoundaryCondition(int matrix_num, unsigned long int index, int dir, int jump);
  void simulation::printAcceptance() const
  {
    if(nAccepts + nRejects > 0)
      std::cout << static_cast<float>(nAccepts)/ static_cast<float>(nAccepts+nRejects) <<std::endl;
    else
      std::cout << "There have been no acceptance or rejections.\n";
  }
const matrix_complex simulation::plaquette(long unsigned index, int dir1, int dir2) const
{
  matrix_complex u1, u2, u3, u4;

  u1 = boundary_condition(dir1, index,dir1,0);
  u2 = boundary_condition(dir2, index,dir1,1);
  u3 = boundary_condition(dir1, index,dir2,1).adjoint();
  u4 = boundary_condition(dir2, index,dir2,0).adjoint();

  return u1*u2*u3*u4;
}
