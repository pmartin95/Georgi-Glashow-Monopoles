#pragma once

#ifndef _LATTICE_
#define __LATTICE__
#include <cstddef>
#include <cstdint>
#include <Eigen/Dense>
#include <random>
#include <vector>
#define DEFAULT_LATTICE_SIZE 4
#define FORALLDIR(i) for(i = 0; i<4; i++)
#define FORALLDIRBUT(i,dir) for(i=0; i<4; i++) if(i != dir)
#define FORALLDIRLESSTHAN(i,j) for( j = 0; j<4; j++) for(i=0; i<j; i++)

typedef Eigen::Matrix2cd matrix_complex;
using namespace std;
class simulation;
//==============================================================================
//LATTICE SITE CLASS============================================================
class lattice_site {
public:
friend class lattice_base;
friend class simulation;
//Constructors------------------------------------------------------------------
lattice_site();
lattice_site(const matrix_complex& link0,const matrix_complex& link1, const matrix_complex& link2, const matrix_complex& link3, const matrix_complex& higg_temp  );
lattice_site(const matrix_complex& link_temp, const matrix_complex& higg_temp  );

lattice_site(const lattice_site& site1); //copy
lattice_site(lattice_site&& site1); //move

//Destructor--------------------------------------------------------------------
~lattice_site(){
}
//Initialization----------------------------------------------------------------
void init_lattice_site();
void init_plattice_site();
void init_lattice_site(std::mt19937_64& g);
void init_plattice_site(std::mt19937_64& g);
//Assignment operators----------------------------------------------------------
lattice_site& operator =(const lattice_site& site1); //copy
lattice_site& operator =( lattice_site&& site1); //move
const matrix_complex& output() const;
const matrix_complex& output(int i) const;
friend std::ostream& operator <<(std::ostream& outputStream,const lattice_site& site1);
matrix_complex link[4];
matrix_complex higgs;
};

//==============================================================================
//LATTICE BASE CLASS============================================================
class lattice_base {
public:
friend class simulation;
friend class Plattice;
//Setup functions
lattice_base();
//Copy Constructor
lattice_base(const lattice_base & L_in);
//Move Constructor
lattice_base(lattice_base&& L_in);
lattice_base(int Nt, int Nx, int Ny, int Nz);
lattice_base(const lattice_site& site1);
//Copy assignment operator
lattice_base& operator =(const lattice_base& L_in);
//Move assignment operator
lattice_base& operator =(lattice_base&& L_in);
~lattice_base();
//Utilty functions
void setLatticeSite(long unsigned int site_index, const matrix_complex& A);
void setLatticeSite(long unsigned int site_index,int dir, const matrix_complex& A);
void infoPrint() const;
long unsigned int coordinateToIndex(int t,int x, int y, int z) const;
long unsigned int coordinateToIndex(const int (&coordinates)[4]) const;
long unsigned int jumpIndex(long unsigned int index, int dir, int jump);
long unsigned int jumpIndex(long unsigned int index,int jump[4]);
void indexToCoordinate(long unsigned int index,int (&coordinates)[4]) const;
void indexToCoordinate(long unsigned int index,int& t, int& x, int& y, int& z) const;
//
const matrix_complex directMatCall(unsigned long int site_index,int matrix_num) const;
int shiftToLattice(int coordinate, int dir) const;
long unsigned shiftToLattice(const long unsigned site_index,const int jump[4]) const;
int incrementCoordinate(int coordinate,int dir) const;//
//private:
int nt,nx,ny,nz;
int ns[4];
long int nsites;
lattice_site * site;
};
//==============================================================================
//LATTICE DERIVED CLASS=========================================================
class lattice : public lattice_base
{
public:
friend class simulation;
friend class Plattice;
//Setup functions
lattice();
lattice(int Nt, int Nx, int Ny, int Nz);
lattice(std::mt19937_64& g);
lattice(int Nt, int Nx, int Ny, int Nz,std::mt19937_64& g);
};
//==============================================================================
//MOMENTUM LATTICE CLASS========================================================
class Plattice : public lattice_base {
public:
friend class simulation;
friend class lattice;
//Setup functions
Plattice();
Plattice(std::mt19937_64& g);
Plattice(std::mt19937_64& g, const lattice& L_in);   //
Plattice(int Nt, int Nx, int Ny, int Nz,std::mt19937_64& g);
};

bool isJump(const int jump[4]);

#endif
