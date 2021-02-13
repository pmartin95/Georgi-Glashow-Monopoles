#pragma once

#ifndef _LATTICE_
#define __LATTICE__

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

class lattice_site {
public:
friend class lattice;
friend class simulation;
lattice_site();
lattice_site(const matrix_complex& link0,const matrix_complex& link1, const matrix_complex& link2, const matrix_complex& link3, const matrix_complex& higg_temp  );
lattice_site(const matrix_complex& link_temp, const matrix_complex& higg_temp  );
lattice_site(const lattice_site& site1);
lattice_site(std::mt19937_64& g);
const matrix_complex& output() const;
const matrix_complex& output(int i) const;
friend std::ostream& operator <<(std::ostream& outputStream,const lattice_site& site1);
lattice_site& operator =(const lattice_site& site1);
//private:
matrix_complex link[4];
matrix_complex higgs;
};

class Plattice_site {
public:
friend class Plattice;
friend class simulation;
Plattice_site();
Plattice_site(const matrix_complex& link0,const matrix_complex& link1, const matrix_complex& link2, const matrix_complex& link3, const matrix_complex& higg_temp  );
Plattice_site(const Plattice_site& site1);
Plattice_site(std::mt19937_64& g);
const matrix_complex& output() const;
const matrix_complex& output(int i) const;
friend std::ostream& operator <<(std::ostream& outputStream,const Plattice_site& site1);
Plattice_site& operator =(const Plattice_site& site1);
//private:
matrix_complex link[4];
matrix_complex higgs;
};

class lattice {
public:
friend class simulation;
friend class Plattice;
//Setup functions
lattice();
lattice(int Nt, int Nx, int Ny, int Nz);
lattice(const lattice& lattice1);
lattice(const lattice_site& site1);
lattice(std::mt19937_64& g);
lattice(int Nt, int Nx, int Ny, int Nz,std::mt19937_64& g);
~lattice();
lattice& operator =(const lattice& L_in);
//Utilty functions
void setLatticeSite(long unsigned int site_index, const matrix_complex& A);
void setLatticeSite(long unsigned int site_index,int dir, const matrix_complex& A);
void infoPrint() const;
long unsigned int coordinateToIndex(int t,int x, int y, int z) const;
long unsigned int coordinateToIndex(int (&coordinates)[4]) const;
long unsigned int jumpIndex(long unsigned int index, int dir, int jump);
long unsigned int jumpIndex(long unsigned int index,int (&jump)[4]);
void indexToCoordinate(long unsigned int index,int (&coordinates)[4]) const;
void indexToCoordinate(long unsigned int index,int& t, int& x, int& y, int& z) const;
//private:
int nt,nx,ny,nz;
int ns[4];
long int nsites;
//std::vector<lattice_site> site;
lattice_site * site;
};

class Plattice {
public:
friend class simulation;
friend class lattice;
//Setup functions
Plattice();
Plattice(const Plattice& Plattice1);
Plattice(const Plattice_site& site1);
Plattice(std::mt19937_64& g);
Plattice(std::mt19937_64& g, const lattice& L_in);   //
Plattice(int Nt, int Nx, int Ny, int Nz,std::mt19937_64& g);
~Plattice();
Plattice& operator =(const Plattice& P_in);
//Utilty functions
void infoPrint() const;
long unsigned int coordinateToIndex(int t,int x, int y, int z) const;
long unsigned int coordinateToIndex(int (&coordinates)[4]) const;
long unsigned int jumpIndex(long unsigned int index, int dir, int jump);
void indexToCoordinate(long unsigned int index,int (&coordinates)[4]) const;
void indexToCoordinate(long unsigned int index,int& t, int& x, int& y, int& z) const;
//private:
int nt,nx,ny,nz;
int ns[4];
long int nsites;
//std::vector<Plattice_site> site;
Plattice_site * site;
};


#endif
