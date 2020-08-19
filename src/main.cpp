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
matrix_complex myExp(const matrix_complex & A, int N);
int main()
{
//   matrix_complex higgsZero, linkIdentity;
//   higgsZero.setZero(); linkIdentity.setIdentity();
//   lattice_site cold_site(linkIdentity, higgsZero);
//   lattice cold_lattice(cold_site);
//
//   // //Initialize Simulation -> lattice
//
//   simulation sim1;
//   simulation sim2(cold_lattice);
//   int stepsArray[3];
//   double HdiffArray1[3] = {0.0d}, HdiffArray2[3] = {0.0d};
//   for(int i = 0; i < 5;i++)
//   {
//     sim1.initializeHMC();
//     sim2.initializeHMC();
//   }
//
//
//   int current_steps = 10;
//   // //Run HMC
//   for(int j=0; j<3;j++)
//   {
//     stepsArray[j] = current_steps;
//     sim1.setupSteps(current_steps);
//     sim2.setupSteps(current_steps);
//     for(int i=0; i < 1;i++)
//     {
//       HdiffArray1[j] += std::abs(sim1.runLeapfrogSimulation());
//       HdiffArray2[j] += std::abs(sim2.runLeapfrogSimulation());
//     }
//     HdiffArray1[j] /= 1.0;
//     HdiffArray2[j] /= 1.0;
//     current_steps *= 10;
//   }
//
// ofstream hot_data, cold_data;
// hot_data.open("hot_data1.txt",std::ios::out);
// cold_data.open("cold_data1.txt",std::ios::out);
// for(int i = 0;i <3; i++)
// {
//   hot_data << stepsArray[i] << " " << HdiffArray1[i] << std::endl;
//   cold_data << stepsArray[i] << " " << HdiffArray2[i] << std::endl;
// }
//
// hot_data.close();
// cold_data.close();
    //std::cout << "Hamiltonian: " << sim1.Hamiltonian() << std::endl;

  //sim1.printAcceptance();

  //Measure observables (average plaquette)
  // std::cout << "<Phi> = " << sim1.averagePhi() << std::endl;
  // std::cout << "<Phi^2> = " << sim1.averagePhi2() << std::endl;
  // std::cout << "Average plaquette = " << sim1.averagePlaquettes() << std::endl;
struct timeval myStartTime;
std::mt19937_64 randTemp(seedGen());

matrix_complex A = uniformSU2Matrix(randTemp);
matrix_complex B;
mytime time1;
std::cout << A << std::endl;
time1.stopwatchStart();
for(int i = 0;i < 11;i++)
  B = A.exp();
std::cout << time1.stopwatchReadSeconds() << std::endl;
std::cout << B << std::endl << std::endl;


time1.stopwatchStart();
for(int i = 0;i < 1;i++)
  B = myExp(A,30);
std::cout << time1.stopwatchReadSeconds() << std::endl;
std::cout << B << std::endl << std::endl;

  return 0;
}


matrix_complex myExp(const matrix_complex & A, int N)
{
  matrix_complex temp, subtotal;
  int factorial;
  int i;
  temp.setIdentity();
  subtotal.setIdentity();
  factorial = 1.0d;
  for(i = 0;i < N;i++)
  {
    temp *=A;
    factorial *= i+1;
    subtotal += temp / static_cast<double>(factorial);
  }
  return subtotal;
}
