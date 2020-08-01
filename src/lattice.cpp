#include <iostream>
#include <complex>
#include "lattice.h"

lattice_site::lattice_site()
{
  for(int i=0;i<4;i++)
    link[i].setIdentity();
  higgs.setZero();
}

lattice_site::lattice_site(const matrix_complex& link0,const matrix_complex& link1, const matrix_complex& link2, const matrix_complex& link3, const matrix_complex& higg_temp  )
{
  link[0] = link0;
  link[1] = link1;
  link[2] = link2;
  link[3] = link3;
  higgs = higg_temp;
}

lattice_site::lattice_site(const lattice_site& site1)
{
  for(int i=0;i<4;i++)
    link[i] = site1.link[i];
  higgs = site1.higgs;
}

lattice_site::lattice_site(URNG& g)
{
  int i;
  FORALLDIR(i)
  link[i] = uniformSU2Matrix(g);
  higgs.setZero();
}

const matrix_complex& lattice_site::output() const
{
  return higgs;
}

const matrix_complex& lattice_site::output(int i) const
{
  return link[i];
}

std::ostream& operator <<(std::ostream& outputStream,const lattice_site& site1)
{
  for(int i=0;i<4;i++)
    outputStream << site1.link[i] << std::endl << std::endl;
  outputStream << site1.higgs;
  return outputStream;
}

lattice_site& lattice_site::operator =(const lattice_site& site1)
{
  if(this == &site1)
  {
    return *this;
  }
  else
  {
    for(int i=0;i<4;i++)
      link[i] = site1.link[i];
    higgs = site1.higgs;
    return *this;
  }
}





lattice::lattice()
{
  nt = DEFAULT_LATTICE_SIZE;
  nx = DEFAULT_LATTICE_SIZE;
  ny = DEFAULT_LATTICE_SIZE;
  nz = DEFAULT_LATTICE_SIZE;
  nsites = nt * nx * ny * nz;
  ns[0] = nt;
  ns[1] = nx;
  ns[2] = ny;
  ns[3] = nz;
  L = new lattice_site[nsites];
  for(long unsigned int i = 0; i < nsites;i++)
  {
    L[i] = lattice_site();
  }
}

lattice::lattice(int Nt, int Nx, int Ny, int Nz)
{
  nt = Nt;
  nx = Nx;
  ny = Ny;
  nz = Nz;

  ns[0] = nt;
  ns[1] = nx;
  ns[2] = ny;
  ns[3] = nz;

  nsites = nt * nx * ny * nz;

  L = new lattice_site[nsites];
  for(long unsigned int i = 0; i < nsites;i++)
  {
    L[i] = lattice_site();
  }
}

lattice::lattice(const lattice& lattice1)
{
  nt = lattice1.nt;
  nx = lattice1.nx;
  ny = lattice1.ny;
  nz = lattice1.nz;
  nsites = nt * nx * ny * nz;
  ns[0] = nt;
  ns[1] = nx;
  ns[2] = ny;
  ns[3] = nz;
  L = new lattice_site[nsites];
  for(long unsigned int i = 0; i < nsites;i++)
  {
    L[i] = lattice1.L[i];
  }
}

lattice::lattice(const lattice_site& site1)
{
  nt = DEFAULT_LATTICE_SIZE;
  nx = DEFAULT_LATTICE_SIZE;
  ny = DEFAULT_LATTICE_SIZE;
  nz = DEFAULT_LATTICE_SIZE;
  nsites = nt * nx * ny * nz;
  ns[0] = nt;
  ns[1] = nx;
  ns[2] = ny;
  ns[3] = nz;
  L = new lattice_site[nsites];
  for(long unsigned int i = 0; i < nsites;i++)
  {
    L[i] = site1;
  }
}

lattice::lattice(URNG& g)
{
  nt = DEFAULT_LATTICE_SIZE;
  nx = DEFAULT_LATTICE_SIZE;
  ny = DEFAULT_LATTICE_SIZE;
  nz = DEFAULT_LATTICE_SIZE;
  nsites = nt * nx * ny * nz;
  ns[0] = nt;
  ns[1] = nx;
  ns[2] = ny;
  ns[3] = nz;
  L = new lattice_site[nsites];
  for(long unsigned int i = 0; i < nsites;i++)
  {
    L[i] = lattice_site(URNG& g);
  }
}


lattice::lattice(int Nt, int Nx, int Ny, int Nz,URNG& g)
{
  nt = Nt;
  nx = Nx;
  ny = Ny;
  nz = Nz;

  ns[0] = nt;
  ns[1] = nx;
  ns[2] = ny;
  ns[3] = nz;

  nsites = nt * nx * ny * nz;

  L = new lattice_site[nsites];
  for(long unsigned int i = 0; i < nsites;i++)
  {
    L[i] = lattice_site(URNG& g);
  }
}


lattice::~lattice()
{
  delete [] L;
}

void lattice::infoPrint() const
{
  std::cout << "This is a " << nt << "x" << nx << "x" << ny << "x" << nz << " lattice.\n";
  std::cout << "It has a total number of " << nsites << " lattice sites.\n";
}

long unsigned int lattice::coordinateToIndex(int t,int x, int y, int z) const
{
  long unsigned int temp;
  temp = y + ny * z;
  temp = x + nx * temp;
  temp = t + nt * temp;
  return temp;
}

long unsigned int lattice::coordinateToIndex(int (&coordinates)[4]) const
{
  int t,x,y,z;
  t = coordinates[0];
  x = coordinates[1];
  y = coordinates[2];
  z = coordinates[3];
  long unsigned int temp;
  temp = y + ny * z;
  temp = x + nx * temp;
  temp = t + nt * temp;
  return temp;
}

long unsigned int lattice::jumpIndex(long unsigned int index, int dir, int jump)
{
  int x[4] = {0};
  indexToCoordinate(index, x);
  x[dir] += jump;
  return coordinateToIndex(x);
}

void lattice::indexToCoordinate(long unsigned int index,int (&coordinates)[4]) const
{
  long unsigned int temp = index;
  coordinates[0] = temp%nt;
  temp = (temp - coordinates[0])/nt;

  coordinates[1] = temp%nx;
  temp = (temp - coordinates[1])/nx;

  coordinates[2] = temp%ny;
  temp = (temp - coordinates[2])/ny;

  coordinates[3] = temp;
}

void lattice::indexToCoordinate(long unsigned int index,int& t, int& x, int& y, int& z) const
{
  long unsigned int temp = index;
  t = temp%nt;
  temp = (temp - t)/nt;

  x = temp%nx;
  temp = (temp - x)/nx;

  y = temp%ny;
  temp = (temp - y)/ny;

  z = static_cast<int>(temp);
}
