#include<iostream>
#include "lattice.h"
#include "rand.h"
#include "sim.h"

int main()
{

  //Initialize Simulation -> lattice

  simulation sim1;
  sim1.printSite(10);
  sim1.printDerivatives(10);
  for(int i=0; i < 10;i++)
    sim1.runLeapfrogSimulation();
  sim1.printSite(10);
    //std::cout << "Hamiltonian: " << sim1.Hamiltonian() << std::endl;
  sim1.printAcceptance();
  std::cout << "<Phi> = " << sim1.averagePhi() << std::endl;
  //Run HMC
  //Measure observables (average plaquette)
  std::cout << "Here.\n";
}
