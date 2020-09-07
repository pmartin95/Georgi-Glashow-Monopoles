#include<iostream>
#include<algorithm>
#include<complex>
#include<math.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <omp.h>
#include "sim.h"
#include "stopwatch.h"

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
  L_temp = Lcopy;
  startMomentum = Plattice(randomGenerator);
  endMomentum = startMomentum;
  P_temp = endMomentum;
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

simulation::~simulation()
{

}

//HMC routines
void simulation::initializeHMC()
{
 for(int i = 0;i< steps;i++)
 {
   leapfrogOneStep();
 }
  std::cout << "\\Delta H: " << georgiGlashowHamiltonian(Lcopy,endMomentum) - georgiGlashowHamiltonian(L,startMomentum) << std::endl;
  //resetMomenta();
  L = Lcopy;
  startMomentum = endMomentum;
}

double simulation::runLeapfrogSimulation()
{
  for(int i = 0;i<steps;i++)
    leapfrogOneStep();

  double Hdiff = georgiGlashowHamiltonian(L, startMomentum) - georgiGlashowHamiltonian(Lcopy, endMomentum);

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
  return Hdiff;
}

void simulation::leapfrogOneStep()
{
  long unsigned int site_index;
  int dir;


  // #pragma omp parallel private(dir) default(none)
  {
    // #pragma omp for
    for(site_index=0; site_index < L.nsites;site_index++)
    {
      //Change the p sign to + instead of -, just to see what happens
      FORALLDIR(dir)
      {
        P_temp.site[site_index].link[dir] = endMomentum.site[site_index].link[dir] - stepSize/2.0 * georgiGlashowActionLinkDerivative(site_index, dir, Lcopy)  ;
      }

      // P_temp.site[site_index].higgs = endMomentum.site[site_index].higgs - stepSize/2.0 * georgiGlashowActionPhiDerivative(site_index, Lcopy);
    }
    // #pragma omp barrier

    // #pragma omp for
    for(site_index=0; site_index < L.nsites;site_index++)
    {
      FORALLDIR(dir)
      L_temp.site[site_index].link[dir] = CayleyHamiltonExp((complex<double>(0.0,1.0) * stepSize * P_temp.site[site_index].link[dir] )) * Lcopy.site[site_index].link[dir] ;
      // L_temp.site[site_index].higgs = Lcopy.site[site_index].higgs + stepSize * P_temp.site[site_index].higgs;
    }
    // #pragma omp barrier

    // if(omp_get_thread_num() == 0)
      Lcopy = L_temp;
    // #pragma omp barrier

    // #pragma omp for
    for(site_index=0; site_index < L.nsites;site_index++)
    {
      FORALLDIR(dir)
      endMomentum.site[site_index].link[dir] = P_temp.site[site_index].link[dir] - stepSize/2.0 * georgiGlashowActionLinkDerivative(site_index, dir, L_temp);
      // endMomentum.site[site_index].higgs = P_temp.site[site_index].higgs - stepSize/2.0 * georgiGlashowActionPhiDerivative(site_index, L_temp);
    }
  }
}

bool simulation::metropolisDecision()
{
  double expResult;
  double randomDecider;
  double H_new, H_old;
  //should return true if to accept configuration, otherwise replace with copy and start over
  H_new = georgiGlashowHamiltonian(Lcopy,endMomentum);
  H_old = georgiGlashowHamiltonian(L,startMomentum);
  std::cout << "H new: " << H_new << " H old: " << H_old << std::endl;
  std::cout << "delta H: " << H_new-H_old << std::endl;
  expResult = std::min(exp(H_old - H_new), 1.0);
  randomDecider = uniformReal(randomGenerator);
  //std::cout << "exp: " << expResult << "  random: " << randomDecider << std::endl;
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
////NOTE FOR PAUL
// I think the issue lies in the phiDerivativePart
//Action functions
double simulation::georgiGlashowLagrangianDensity(long unsigned int site_index) const
{
  //Declaration
  matrix_complex mainPhi, mainLink[4];
  double phiSquareTrace;
  double phiDerivativePart = 0.0, plaquettePart = 0.0, miscPhiPart = 0.0;
  int jumpNone[4] = {0};
  int dir,i,j;
  //Initialization
  mainPhi = matCall(L,5,site_index,jumpNone);
  FORALLDIR(dir)
    mainLink[dir] = matCall(L,dir,site_index,jumpNone);
  phiSquareTrace = (mainPhi * mainPhi).trace().real();

  //Phi "derivative"
  phiDerivativePart = 4.0 * phiSquareTrace;
  FORALLDIR(dir)
  {
    int jumpTemp[4] = {0};
    jumpTemp[dir]++;
    phiDerivativePart -= (mainPhi * mainLink[dir] * matCall(L,dir,site_index,jumpTemp) * mainLink[dir].adjoint() ).trace().real();
  }
  phiDerivativePart *= 2.0;
  //plaquette
  FORALLDIRLESSTHAN(i,j)
  {
    int muJump[4] = {0}, nuJump[4] = {0};
    muJump[i]++; nuJump[j]++;
    plaquettePart += 2.0 - ( mainLink[i] * matCall(L,i,site_index,muJump) * matCall(L,j,site_index,nuJump).adjoint() * mainLink[j].adjoint() ).trace().real();
  }
  plaquettePart *= 2.0/(g*g);
  //other phi terms
  miscPhiPart = m2 * phiSquareTrace + lambda * phiSquareTrace *phiSquareTrace;
  phiDerivativePart = 0.0; ////
  return phiDerivativePart + plaquettePart + miscPhiPart;
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
    for(long unsigned int i = 0; i<L_in.nsites;i++)
    {
      local_total += georgiGlashowLagrangianDensity(i,L_in);
    }
    #pragma omp critical
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
    matrix_complex local_momenta;
    local_momenta.setZero();
    int i;
    #pragma omp for
    for(site_index=0;site_index < P_in.nsites; site_index++)
    {
      FORALLDIR(i)
        local_momenta+= P_in.site[site_index].link[i] *  P_in.site[site_index].link[i];
      local_momenta+=   P_in.site[site_index].higgs *  P_in.site[site_index].higgs;
    }
    #pragma omp critical
    momenta_matrix += local_momenta;
  }
  momenta_total = momenta_matrix.trace().real();
  return momenta_total/4.0;
}


double simulation::georgiGlashowHamiltonian(const lattice& L_in, const Plattice& P_in) const
{
  matrix_complex momenta_matrix;
  double field_total, kinetic_total;
  long unsigned int site_index;
  int i;

  field_total = georgiGlashowAction(L_in);
    std::cout << "Field total: " << field_total << std::endl;

  kinetic_total = kineticTerm(P_in);
    std::cout << "Momenta total: " << kinetic_total << std::endl;
  return field_total + kinetic_total;
}

double simulation::Hamiltonian() const
{
  return georgiGlashowHamiltonian(L,startMomentum) /L.nsites ;
}

const matrix_complex simulation::georgiGlashowActionLinkDerivative(long unsigned int site_index, int dir, const lattice& L_in) const //
{
  matrix_complex topStaple, bottomStaple;
  matrix_complex temp1,temp2;
  matrix_complex Ulink;
  matrix_complex subtotal1,subtotal2;
  matrix_complex identityMatrix;

  int jump1[4] = {0};
  int jumpNone[4] = {0,0,0,0};
  int temp_dir;
  identityMatrix.setIdentity();

  Ulink = matCall(L_in,dir,site_index,jumpNone);
  std::cout << "Ulink: " << Ulink << std::endl;
  topStaple.setZero(); bottomStaple.setZero();
  temp1.setZero(); temp2.setZero();
  subtotal1.setZero(); subtotal2.setZero();
  jump1[dir]++;

  FORALLDIRBUT(temp_dir,dir)
  {
      int jump2[4] = {0}, jump3[4] = {0}, jump4[4] = {0};
      jump2[temp_dir]++; jump3[temp_dir]--;
      jump4[dir]++; jump4[temp_dir]--;

      topStaple.noalias() += matCall(L_in, site_index,temp_dir,jump1) * matCall(L_in, site_index,dir,jump2).adjoint() * matCall(L_in, site_index,temp_dir,jumpNone).adjoint();
      bottomStaple.noalias() += matCall(L_in, site_index,temp_dir,jump4).adjoint() * matCall(L_in, site_index,dir,jump3).adjoint() * matCall(L_in, site_index,temp_dir,jump3);
  }

  subtotal1.noalias() += complex<double>(0.0,1.0/(g*g)) * (  Ulink*(topStaple + bottomStaple)  - (topStaple + bottomStaple).adjoint()  * Ulink.adjoint()   );
  subtotal1 = subtotal1 - subtotal1.trace()/2.0 * identityMatrix;
  temp1 = Ulink * matCall(L_in,4,site_index,jump1) * matCall(L_in,dir,site_index,jumpNone).adjoint()  ;
  //std::cout << "Link force, temp1: " << temp1 << std::endl;
  temp2 = matCall(L_in,4,site_index,jumpNone);
  //std::cout << "Link force, temp2: " << temp1 << std::endl;
  subtotal2.noalias() += (temp1*temp2 - temp2*temp1) * complex<double>(0.0,2.0);
  //std::cout << "Subtotal1 trace: " << subtotal1.trace() << std::endl;
  //std::cout << "Subtotal2 trace: " << subtotal2.trace() << std::endl;

  return subtotal1 + subtotal2 ;// best so far
}

const matrix_complex simulation::georgiGlashowActionPhiDerivative(long unsigned int site_index, const lattice& L_in) const
{
    matrix_complex temp1, temp2, temp3,Iden,temp;
    int temp_dir;
    int jumpNone[4] = {0,0,0,0};
    double traced_part;
    temp2.setZero();
    temp3.setZero();
    Iden.setIdentity();

    // FORALLDIR(temp_dir)
    // {
    //   int jump1[4] = {0}, jump2[4] ={0};
    //   jump1[temp_dir]++;jump2[temp_dir]--;
    //   temp2.noalias() += matCall(L_in,temp_dir,site_index,jumpNone) * matCall(L_in,4,site_index,jump1) * matCall(L_in,temp_dir,site_index,jumpNone).adjoint();
    //   temp3.noalias() += matCall(L_in,temp_dir,site_index,jump2).adjoint() *  matCall(L_in,4,site_index,jump2)* matCall(L_in,temp_dir,site_index,jump2);
    // }
    temp1 = matCall(L_in,4,site_index,jumpNone);
    temp1.noalias() = (8.0 * lambda *temp1*(temp1*temp1).trace() + (4.0*m2*+ 32.0)*temp1  );

    //std::cout << "Link force, temp1: " << temp1 << std::endl;
    temp = 4.0*(temp2 + temp3) - temp1 ;
    temp = temp - temp.trace() / 2.0 * Iden;
    return temp;
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
    for(long unsigned int i = 0; i< L.nsites;i++)
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
    for(long unsigned int i = 0; i< L.nsites;i++)
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
    for(long unsigned int i = 0; i< L.nsites;i++)
    {
        local_subtotal+= L.site[i].higgs * L.site[i].higgs;
    }
    #pragma omp critical
    subtotal+=local_subtotal;
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
  startMomentum = Plattice(randomGenerator);
  endMomentum = startMomentum;
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
void simulation::printAcceptance() const
{
  if(nAccepts + nRejects > 0)
  {
    std::cout << static_cast<float>(nAccepts)/ static_cast<float>(nAccepts+nRejects) * 100 << "%\n";
    std::cout << nAccepts << " accepts, " << nRejects << " rejects.\n";
  }
  else
    std::cout << "There have been no acceptances or rejections.\n";
}

void simulation::printSite(long unsigned int site_index) const
{
  std::cout << L.site[site_index].link[0] <<std::endl;
  std::cout << L.site[site_index].link[1] <<std::endl;
  std::cout << L.site[site_index].link[2] <<std::endl;
  std::cout << L.site[site_index].link[3] <<std::endl;
  std::cout << L.site[site_index].higgs <<std::endl <<std::endl;
}

void simulation::printDerivatives(long unsigned int site_index) const
{
  std::cout << "Calling print derivative routine\n\n.";
  int dir;
  FORALLDIR(dir)
    std::cout << georgiGlashowActionLinkDerivative(site_index,dir, L)<< std::endl;
  std::cout << georgiGlashowActionPhiDerivative(site_index, L) << std::endl << std::endl;
  std::cout << "Exiting print derivative routine\n\n.";
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


matrix_complex CayleyHamiltonExp(const matrix_complex & A)
{

  complex<double> M = std::sqrt( A(0,1) *A(1,0) - A(0,0) *A(1,1)   );
  matrix_complex iden;
  iden.setIdentity();

  return std::cosh(M) * iden + std::sinh(M)/M *A   ;
}
