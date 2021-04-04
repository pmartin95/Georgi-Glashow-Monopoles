#include <iostream>
#include <complex>
#include <cstdlib>
#include <Eigen/Dense>
#include "lattice.h"
#include "rand.h"

Plattice_site::Plattice_site()
{
        int i;
        FORALLDIR(i)
        link[i].setZero();
        higgs.setZero();
}

Plattice_site::Plattice_site(std::mt19937_64& g)
{
        int i;
        FORALLDIR(i)
        link[i] = normalHermitianMatrix(g);
        higgs = normalHermitianMatrix(g);
}
//copy constructor
Plattice_site::Plattice_site(const Plattice_site& site1)
{
        for(int i=0; i<4; i++)
                link[i] = site1.link[i];
        higgs = site1.higgs;
}

//move constructor
Plattice_site::Plattice_site(Plattice_site&& site1)
{
        int i;
        FORALLDIR(i)
        {
                link[i] = site1.link[i];
        }
        higgs = site1.higgs;
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

Plattice_site& Plattice_site::operator =(Plattice_site&& site1)
{
        if (this != &site1)
        {
                int i;
                FORALLDIR(i)
                {
                        link[i] = site1.link[i];
                }
                higgs = site1.higgs;
        }
        return *this;
}
