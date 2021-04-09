
#include <iostream>
#include <complex>
#include <cstdlib>
#include <Eigen/Dense>
#include "lattice.h"
#include "rand.h"


//////// Here below needs ot be looked over
// Need to redo functions for P lattice sites



Plattice::Plattice()
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
        site = new Plattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = Plattice_site();
        }
}


Plattice::Plattice(std::mt19937_64& g)
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
        site = new Plattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = Plattice_site(g);
        }
}

Plattice::Plattice(std::mt19937_64& g, const lattice& L_in)
{
        nt = L_in.nt;
        nx = L_in.nx;
        ny = L_in.ny;
        nz = L_in.nz;
        nsites = nt * nx * ny * nz;
        ns[0] = nt;
        ns[1] = nx;
        ns[2] = ny;
        ns[3] = nz;
        site = new Plattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = Plattice_site(g);
        }
}

Plattice::Plattice(int Nt, int Nx, int Ny, int Nz,std::mt19937_64& g)
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

        site = new Plattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = Plattice_site(g);
        }
}
//Copy Constructor
Plattice::Plattice(const Plattice& P_in)
{
        nt = P_in.nt;
        nx = P_in.nx;
        ny = P_in.ny;
        nz = P_in.nz;
        nsites = nt * nx * ny * nz;
        ns[0] = nt;
        ns[1] = nx;
        ns[2] = ny;
        ns[3] = nz;
        site = new Plattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = P_in.site[i];
        }
}
//Move Constructor
Plattice::Plattice(Plattice&& P_in) : site(nullptr)
{
        nt = P_in.nt;
        nx = P_in.nx;
        ny = P_in.ny;
        nz = P_in.nz;
        nsites = nt * nx * ny * nz;
        ns[0] = nt;
        ns[1] = nx;
        ns[2] = ny;
        ns[3] = nz;
        site = P_in.site;
        P_in.site = nullptr;
}
//Copy assignment
Plattice& Plattice::operator =(const Plattice& P_in)
{
        if(this != &P_in)
        {
                if( this->nsites == P_in.nsites &&  this->nt == P_in.nt &&  this->nx == P_in.nx &&  this->ny == P_in.ny &&  this->nz == P_in.nz )
                {
                        for(long unsigned int i = 0; i < nsites; i++)
                        {
                                site[i] = P_in.site[i];
                        }
                }
                else
                {
                        nt = P_in.nt;
                        nx = P_in.nx;
                        ny = P_in.ny;
                        nz = P_in.nz;

                        ns[0] = nt;
                        ns[1] = nx;
                        ns[2] = ny;
                        ns[3] = nz;

                        nsites = nt * nx * ny * nz;
                        delete [] site;
                        site = new Plattice_site[nsites];
                        for(long unsigned int i = 0; i < nsites; i++)
                        {
                                site[i] = P_in.site[i];
                        }
                }
        }
        return *this;
}

Plattice& Plattice::operator =(Plattice&& P_in)
{
        if(this != &P_in)
        {
                delete[] site;
                nt = P_in.nt;
                nx = P_in.nx;
                ny = P_in.ny;
                nz = P_in.nz;

                ns[0] = nt;
                ns[1] = nx;
                ns[2] = ny;
                ns[3] = nz;

                nsites = nt * nx * ny * nz;
                site = P_in.site;
                P_in.site = nullptr;
        }
        return *this;
}
Plattice::~Plattice()
{
        delete [] site;
}


void Plattice::infoPrint() const
{
        std::cout << "This is a " << nt << "x" << nx << "x" << ny << "x" << nz << " lattice.\n";
        std::cout << "It has a total number of " << nsites << " lattice sites.\n";
}

long unsigned int Plattice::coordinateToIndex(int t,int x, int y, int z) const
{
        long unsigned int temp;
        temp = y + ny * z;
        temp = x + nx * temp;
        temp = t + nt * temp;
        return temp;
}

long unsigned int Plattice::coordinateToIndex(int coordinates[4]) const
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

long unsigned int Plattice::jumpIndex(long unsigned int index, int dir, int jump)
{
        int x[4] = {0};
        indexToCoordinate(index, x);
        x[dir] += jump;
        return coordinateToIndex(x);
}

void Plattice::indexToCoordinate(long unsigned int index,int (&coordinates)[4]) const
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

void Plattice::indexToCoordinate(long unsigned int index,int& t, int& x, int& y, int& z) const
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
