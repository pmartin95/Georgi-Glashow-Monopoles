#include <cstdlib>
#include <sys/time.h>
#include<iostream>
#include<fstream>
#include<complex>
#include<cstdlib>
#include "lattice.h"
#include "rand.h"
#include "sim.h"
#include "stopwatch.h"
#include <unsupported/Eigen/MatrixFunctions>

int main()
{
  matrix_complex higgsZero, linkIdentity;
  higgsZero.setZero(); linkIdentity.setIdentity();
  lattice_site cold_site(linkIdentity, higgsZero);
  lattice cold_lattice(cold_site);
//
//   // //Initialize Simulation -> lattice
//
  simulation sim1;
  simulation sim2(cold_lattice);

  int stepsArray[3];
  double HdiffArray1[3] = {0.0d}, HdiffArray2[3] = {0.0d};


  std::cout << "force * force.adjoint()= " << sim1.georgiGlashowActionLinkDerivative(10,1,sim1.L) - sim1.georgiGlashowActionLinkDerivative(10,1,sim1.L).adjoint() << std::endl;
  std::cout << "force * force.adjoint()= " << sim1.georgiGlashowActionPhiDerivative(10,sim1.L) - sim1.georgiGlashowActionPhiDerivative(10,sim1.L).adjoint() << std::endl;
  sim1.setupSteps(10);
  sim2.setupSteps(10);
  for(int i = 0; i < 30;i++)
  {
    mytime time1;
    time1.stopwatchStart();
    sim1.initializeHMC();
    std::cout << "Initialization step: " << i << " of 30. Sim1 time: " << time1.stopwatchReadSeconds() << std::endl << std::endl;
    time1.stopwatchStart();
    sim2.initializeHMC();
    std::cout << "Initialization step: " << i << " of 30. Sim2 time: " << time1.stopwatchReadSeconds() << std::endl << std::endl;

  }
//
//
  int current_steps = 10;
  // //Run HMC
  for(int j=0; j<3;j++)
  {
    stepsArray[j] = current_steps;
    sim1.setupSteps(current_steps);
    sim2.setupSteps(current_steps);
    for(int i=0; i < 30;i++)
    {
      HdiffArray1[j] += std::abs(sim1.runLeapfrogSimulation());
      HdiffArray2[j] += std::abs(sim2.runLeapfrogSimulation());
    }
    HdiffArray1[j] /= 30.00;
    HdiffArray2[j] /= 30.00;
    current_steps *= 10;
  }

ofstream hot_data, cold_data;
hot_data.open("hot_data1.txt",std::ios::out);
cold_data.open("cold_data1.txt",std::ios::out);
for(int i = 0;i <3; i++)
{
  hot_data << stepsArray[i] << " " << HdiffArray1[i] << std::endl;
  cold_data << stepsArray[i] << " " << HdiffArray2[i] << std::endl;
}

hot_data.close();
cold_data.close();
    std::cout << "Hamiltonian: " << sim1.Hamiltonian() << std::endl;

  sim1.printAcceptance();
  sim1.printAcceptance();

  ////Measure observables (average plaquette)
  std::cout << "<Phi> = " << sim1.averagePhi() << std::endl;
  std::cout << "<Phi^2> = " << sim1.averagePhi2() << std::endl;
  std::cout << "Average plaquette = " << sim1.averagePlaquettes() << std::endl;

// simulation sim1;
//
// for(int i = 10;i < 100;i*=10 )
// {
//   sim1.setupSteps(i);
//   std::cout << "i= " << i << std::endl;
//   sim1.runLeapfrogSimulation();
//   std::cout  << std::endl << std::endl;
// }


  return 0;
}
