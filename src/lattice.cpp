#include <iostream>
#include <complex>
#include <cstdlib>
#include <Eigen/Dense>
#include "lattice.h"
#include "rand.h"

lattice_site::lattice_site()
{
        for(int i=0; i<4; i++)
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
lattice_site::lattice_site(const matrix_complex& link_temp, const matrix_complex& higg_temp  )
{
        link[0] = link_temp;
        link[1] = link_temp;
        link[2] = link_temp;
        link[3] = link_temp;
        higgs = higg_temp;
}


lattice_site::lattice_site(const lattice_site& site1)
{
        for(int i=0; i<4; i++)
                link[i] = site1.link[i];
        higgs = site1.higgs;
}

lattice_site::lattice_site(std::mt19937_64& g)
{
        int i;
        FORALLDIR(i)
        link[i] = uniformSU2Matrix(g);
        higgs = 0.005* normalHermitianMatrix(g);
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
        for(int i=0; i<4; i++)
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
                for(int i=0; i<4; i++)
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
        site = new lattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = lattice_site();
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

        site = new lattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = lattice_site();
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
        site = new lattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = lattice1.site[i];
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

lattice& lattice::operator =(const lattice& L_in)
{
        if(this == &L_in)
        {
                return *this;
        }
        else if( this->nsites == L_in.nsites &&  this->nt == L_in.nt &&  this->nx == L_in.nx &&  this->ny == L_in.ny &&  this->nz == L_in.nz )
        {
                for(long unsigned int i = 0; i < nsites; i++)
                {
                        site[i] = L_in.site[i];
                }
                return *this;
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
                return *this;
        }
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

//////// Here below needs ot be looked over
// Need to redo functions for P lattice sites
Plattice_site::Plattice_site()
{
        int i;
        FORALLDIR(i)
        link[i].setZero();
        higgs.setZero();
}

Plattice_site::Plattice_site(const matrix_complex& link0,const matrix_complex& link1, const matrix_complex& link2, const matrix_complex& link3, const matrix_complex& higg_temp  )
{
        link[0] = link0;
        link[1] = link1;
        link[2] = link2;
        link[3] = link3;
        higgs = higg_temp;
}

Plattice_site::Plattice_site(const Plattice_site& site1)
{
        for(int i=0; i<4; i++)
                link[i] = site1.link[i];
        higgs = site1.higgs;
}

Plattice_site::Plattice_site(std::mt19937_64& g)
{
        int i;
        FORALLDIR(i)
        link[i] = normalHermitianMatrix(g);
        higgs = normalHermitianMatrix(g);
}

const matrix_complex& Plattice_site::output() const
{
        return higgs;
}

const matrix_complex& Plattice_site::output(int i) const
{
        return link[i];
}

std::ostream& operator <<(std::ostream& outputStream,const Plattice_site& site1)
{
        for(int i=0; i<4; i++)
                outputStream << site1.link[i] << std::endl << std::endl;
        outputStream << site1.higgs;
        return outputStream;
}

Plattice_site& Plattice_site::operator =(const Plattice_site& site1)
{
        if(this == &site1)
        {
                return *this;
        }
        else
        {
                for(int i=0; i<4; i++)
                        link[i] = site1.link[i];
                higgs = site1.higgs;
                return *this;
        }
}


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

Plattice::Plattice(const Plattice& Plattice1)
{
        nt = Plattice1.nt;
        nx = Plattice1.nx;
        ny = Plattice1.ny;
        nz = Plattice1.nz;
        nsites = nt * nx * ny * nz;
        ns[0] = nt;
        ns[1] = nx;
        ns[2] = ny;
        ns[3] = nz;
        site = new Plattice_site[nsites];
        for(long unsigned int i = 0; i < nsites; i++)
        {
                site[i] = Plattice1.site[i];
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


Plattice::~Plattice()
{
        delete [] site;
}

Plattice& Plattice::operator =(const Plattice& P_in)
{
        if(this == &P_in)
        {
                return *this;
        }
        else if( this->nsites == P_in.nsites &&  this->nt == P_in.nt &&  this->nx == P_in.nx &&  this->ny == P_in.ny &&  this->nz == P_in.nz )
        {
                for(long unsigned int i = 0; i < nsites; i++)
                {
                        site[i] = P_in.site[i];
                }
                return *this;
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
                return *this;
        }
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

long unsigned int Plattice::coordinateToIndex(int (&coordinates)[4]) const
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
