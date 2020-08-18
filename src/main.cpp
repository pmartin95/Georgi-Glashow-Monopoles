#include<iostream>
#include<complex>
#include "lattice.h"
#include "rand.h"
#include "sim.h"
#include <unsupported/Eigen/MatrixFunctions>

int main()
{

//set higgs to zero, link fields to diag{e^{i\phi},e^{-i\phi}}
// matrix_complex A,B;
// double phi= 3.141592/4.0;
// complex<double> temp;
// temp = std::exp( complex<double>(0.0d,phi)  );
// A(0,0) = temp;
// A(1,1) = std::conj(temp);
// B.setZero();
//
// lattice_site my_site(A,B);
// lattice my_L(my_site);
// simulation my_sim(my_L);
//
// std::cout << "Higgs force: " << my_sim.georgiGlashowActionPhiDerivative(1,my_L);
// std::cout << "Link force: " << my_sim.georgiGlashowActionLinkDerivative(1,1,my_L);

  // //Initialize Simulation -> lattice
  simulation sim1;
  sim1.initializeHMC();
  //
  // //Run HMC
  for(int i=0; i < 10;i++)
    sim1.runLeapfrogSimulation();
    //std::cout << "Hamiltonian: " << sim1.Hamiltonian() << std::endl;
  sim1.printAcceptance();

  //Measure observables (average plaquette)
  std::cout << "<Phi> = " << sim1.averagePhi() << std::endl;
  std::cout << "<Phi^2> = " << sim1.averagePhi2() << std::endl;
  std::cout << "Average plaquette = " << sim1.averagePlaquettes() << std::endl;
}
