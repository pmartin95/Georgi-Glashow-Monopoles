#include <random>
#include "rand.h"


unsigned seedGen()
{
  typedef std::chrono::high_resolution_clock myclock;
  myclock::time_point beginning = myclock::now();
  myclock::duration d = myclock::now() - beginning;
  unsigned seed2 = d.count();
  return seed2;
}


const matrix_complex normalHermitianMatrix( URNG& g ) //Also traceless
{
  double a,b,c;
  std::normal_distribution<double> distribution(0.0,1.0);
  matrix_complex herm;
  a = distribution(g);
  b = distribution(g);
  c = distribution(g);
  herm(0,0) = complex(c);
  herm(1,1) = complex(-c);
  herm(0,1) = complex(a,-b);
  herm(1,0) = complex(a,b);
  return herm;
}

const matrix_complex uniformSU2Matrix( URNG& g ) // Identity plus an anti-Hermitian matrix
{
  double a[4], total;
  matrix_complex su;
  total = 0;
  for(int i=0;i<4;i++)
  {
    a[i] = uniformReal(g,-1.0,1.0);
    total += a[i] * a[i];
  }
  total = sqrt(total);
  for(int i=0;i<4;i++)
    a[i] /= total;
  su(0,0) = complex(a[0],a[1]);
  su(1,1) = complex(a[0],-a[1]);
  su(0,1) = complex(a[2],a[3]);
  su(1,0) = complex(-a[2],a[3]);
  return su;
}

const double uniformReal(URNG& g)
{
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  return distribution(g);
}

const double uniformReal(URNG& g, double low, double high)
{
  std::uniform_real_distribution<double> distribution(low,high);
  return distribution(g);
}
