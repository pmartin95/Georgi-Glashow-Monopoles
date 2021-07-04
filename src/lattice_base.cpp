#include <iostream>
#include <complex>
#include <cstdlib>
#include <Eigen/Dense>
#include "lattice.h"
#include "rand.h"


lattice_base::lattice_base()
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
}

//Copy Constructor
lattice_base::lattice_base(const lattice_base & L_in)
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
lattice_base::lattice_base(lattice_base&& L_in) : site(nullptr)
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

lattice_base::lattice_base(int Nt, int Nx, int Ny, int Nz)
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
}


lattice_base::lattice_base(const lattice_site& site1)
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
                site[i] = site1;
}

lattice_base::~lattice_base()
{
        delete [] site;
}
//Copy assignment
lattice_base& lattice_base::operator =(const lattice_base& L_in)
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

lattice_base& lattice_base::operator =(lattice_base&& L_in)
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

void lattice_base::setLatticeSite(long unsigned int site_index, const matrix_complex& A)
{
        site[site_index].higgs = A;
}
void lattice_base::setLatticeSite(long unsigned int site_index,int dir, const matrix_complex& A)
{
        site[site_index].link[dir] = A;
}
void lattice_base::infoPrint() const
{
        std::cout << "This is a " << nt << "x" << nx << "x" << ny << "x" << nz << " lattice.\n";
        std::cout << "It has a total number of " << nsites << " lattice sites.\n";
}

long unsigned int lattice_base::coordinateToIndex(int t,int x, int y, int z) const
{
        long unsigned int temp;
        temp = y + ny * z;
        temp = x + nx * temp;
        temp = t + nt * temp;
        return temp;
}

long unsigned int lattice_base::coordinateToIndex(const int (&coordinates)[4]) const
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

long unsigned int lattice_base::jumpIndex(long unsigned int index, int dir, int jump)
{
        int x[4] = {0};
        indexToCoordinate(index, x);
        x[dir] += jump;
        return coordinateToIndex(x);
}

long unsigned int lattice_base::jumpIndex(long unsigned int index,int jump[4])
{
        int x[4] = {0},dir;
        indexToCoordinate(index, x);
        FORALLDIR(dir)
        x[dir] += jump[dir];
        return coordinateToIndex(x);
}

void lattice_base::indexToCoordinate(long unsigned int index,int (&coordinates)[4]) const
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

void lattice_base::indexToCoordinate(long unsigned int index,int& t, int& x, int& y, int& z) const
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


const matrix_complex lattice_base::directMatCall(unsigned long int site_index,int matrix_num) const
{

        if(matrix_num<4)
                return site[site_index].link[matrix_num];
        else
                return site[site_index].higgs;
}
//The following function shifts a coordinate to be on lattice via adding and
//reducing based on lattice size
int lattice_base::shiftToLattice(int coordinate, int dir) const
{
        int temp_coordinate = coordinate;
        while(temp_coordinate < 0.0)
                temp_coordinate += ns[dir];
        return (temp_coordinate)%(ns[dir]);
}
long unsigned lattice_base::shiftToLattice(const long unsigned site_index,const int jump[4])const
{
        int x[4] ={0}, dir;
        indexToCoordinate(site_index, x);
        FORALLDIR(dir)
        {
                x[dir] += jump[dir];
                x[dir] = shiftToLattice(x[dir],dir);
        }
        return coordinateToIndex(x);
}
int lattice_base::incrementCoordinate(int coordinate,int dir) const
{
        if(coordinate > ns[dir] )
                return coordinate - ns[dir];
        else if(coordinate <  0)
                return coordinate + ns[dir];
        else
        {
                std::cout << "ERROR: Indices not corresponding" << std::endl;
                std::exit(1);
        }
}
//Checks to see if there is a jump
bool isJump(const int jump[4])
{
        int i;
        bool jumpp;
        FORALLDIR(i)
        {
                jumpp = (jump[i] != 0);
                if(jumpp)
                        break;
        }
        return jumpp;
}
