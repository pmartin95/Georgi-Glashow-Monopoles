
#include <iostream>
#include <complex>
#include <cstdlib>
#include <Eigen/Dense>
#include "lattice.h"
#include "rand.h"


//////// Here below needs ot be looked over
// Need to redo functions for P lattice sites



Plattice::Plattice() : lattice_base()
{
        for(long unsigned int i = 0; i < nsites; i++)
                site[i].init_plattice_site();
}


Plattice::Plattice(std::mt19937_64& g) : lattice_base()
{
        for(long unsigned int i = 0; i < nsites; i++)
                site[i].init_plattice_site(g);

}

Plattice::Plattice(std::mt19937_64& g, const lattice& L_in) : lattice_base(L_in)
{
        for(long unsigned int i = 0; i < nsites; i++)
                site[i].init_plattice_site(g);
}

Plattice::Plattice(int Nt, int Nx, int Ny, int Nz,std::mt19937_64& g) : lattice_base(Nt,Nx,Ny,Nz)
{

        for(long unsigned int i = 0; i < nsites; i++)
                site[i].init_plattice_site(g);
}
