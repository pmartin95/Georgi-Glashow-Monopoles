#include <iostream>
#include <complex>
#include <cstdlib>
#include <Eigen/Dense>
#include "lattice.h"
#include "rand.h"


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
        site = new lattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = lattice_site();
        }
}

//Copy Constructor
lattice::lattice(const lattice& L_in)
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
        site = new lattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = L_in.site[i];
        }
}
//Move Constructor
lattice::lattice(lattice&& L_in) : site(nullptr)
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
        site = L_in.site;
        L_in.site = nullptr;
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

        site = new lattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = lattice_site();
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
        site = new lattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = site1;
        }
}

lattice::lattice(std::mt19937_64& g)
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
        site = new lattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = lattice_site(g);
        }
}


lattice::lattice(int Nt, int Nx, int Ny, int Nz,std::mt19937_64& g)
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

        site = new lattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = lattice_site(g);
        }
}


lattice::~lattice()
{
        delete [] site;
}
//Copy assignment
lattice& lattice::operator =(const lattice& L_in)
{
        //std::cout << "Copying lattice...\n";
        if(this != &L_in)
        {
                if( this->nsites == L_in.nsites &&  this->nt == L_in.nt &&  this->nx == L_in.nx &&  this->ny == L_in.ny &&  this->nz == L_in.nz )
                {
                        for(long unsigned int i = 0; i < nsites; i++)
                        {
                                site[i] = L_in.site[i];
                        }
                }
                else
                {
                        nt = L_in.nt;
                        nx = L_in.nx;
                        ny = L_in.ny;
                        nz = L_in.nz;

                        ns[0] = nt;
                        ns[1] = nx;
                        ns[2] = ny;
                        ns[3] = nz;

                        nsites = nt * nx * ny * nz;
                        delete [] site;
                        site = new lattice_site[nsites];
                        for(long unsigned int i = 0; i < nsites; i++)
                        {
                                site[i] = L_in.site[i];
                        }
                }
        }
        return *this;
}

lattice& lattice::operator =(lattice&& L_in)
{
        //std::cout << "Moving lattice...\n";
        if(this != &L_in)
        {
                delete[] site;
                nt = L_in.nt;
                nx = L_in.nx;
                ny = L_in.ny;
                nz = L_in.nz;

                ns[0] = nt;
                ns[1] = nx;
                ns[2] = ny;
                ns[3] = nz;

                nsites = nt * nx * ny * nz;
                site = L_in.site;
                L_in.site = nullptr;
        }
        return *this;
}

void lattice::setLatticeSite(long unsigned int site_index, const matrix_complex& A)
{
        site[site_index].higgs = A;
}
void lattice::setLatticeSite(long unsigned int site_index,int dir, const matrix_complex& A)
{
        site[site_index].link[dir] = A;
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

long unsigned int lattice::jumpIndex(long unsigned int index,int (&jump)[4])
{
        int x[4] = {0},dir;
        indexToCoordinate(index, x);
        FORALLDIR(dir)
        x[dir] += jump[dir];
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
