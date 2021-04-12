#include <iostream>
#include <complex>
#include <cstdlib>
#include <Eigen/Dense>
#include "lattice.h"
#include "rand.h"


lattice::lattice() : lattice_base()
{
        for(long unsigned int i = 0; i < nsites; i++)
                site[i].init_lattice_site();
}

lattice::lattice(int Nt, int Nx, int Ny, int Nz) : lattice_base(Nt, Nx, Ny, Nz)
{
        for(long unsigned int i = 0; i < nsites; i++)
                site[i] = lattice_site();
}

lattice::lattice(std::mt19937_64& g) : lattice_base()
{

        for(long unsigned int i = 0; i < nsites; i++)
                site[i].init_lattice_site(g);
}

lattice::lattice(int Nt, int Nx, int Ny, int Nz,std::mt19937_64& g) : lattice_base(Nt,Nx,Ny,Nz)
{
        for(long unsigned int i = 0; i < nsites; i++)
                site[i].init_lattice_site(g);
}
