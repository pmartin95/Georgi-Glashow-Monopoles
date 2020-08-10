#include<iostream>
#include "lattice.h"
#include "rand.h"
#include "sim.h"

int main()
{

  // Plattice startMomentum;
  // startMomentum = Plattice();
  // std::cout << "Setup Plattice\n";
  // Plattice endMomentum;
  // endMomentum = startMomentum;
  // std::cout << "Copied Plattice\n";

  int i;
  //Initialize Simulation -> lattice
  simulation sim1;
  for(i = 0; i< 50;i++)
    sim1.runLeapfrogSimulation();
    //std::cout << "Hamiltonian: " << sim1.Hamiltonian() << std::endl;
  sim1.printAcceptance();
  std::cout << "<Phi> = " << sim1.averagePhi() << std::endl;
  //Run HMC
  //Measure observables (average plaquette)
  std::cout << "Here.\n";
}
