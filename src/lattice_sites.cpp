#include <iostream>
#include <complex>
#include <cstdlib>
//#include <cstddef>
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
lattice_site::lattice_site(std::mt19937_64& g)
{
        int i;
        FORALLDIR(i)
        link[i] = uniformSU2Matrix(g);
        higgs = 0.005* normalHermitianMatrix(g);
}

//copy constructor
lattice_site::lattice_site(const lattice_site& site1)
{
        for(int i=0; i<4; i++)
                link[i] = site1.link[i];
        higgs = site1.higgs;
}

//move constructor
lattice_site::lattice_site(lattice_site&& site1)
{
        int i;
        FORALLDIR(i)
        {
                link[i] = site1.link[i];
        }
        higgs = site1.higgs;
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

//Copy assignment
lattice_site& lattice_site::operator =(const lattice_site& site1)
{
        if(this != &site1)
        {
                for(int i=0; i<4; i++)
                        link[i] = site1.link[i];
                higgs = site1.higgs;
        }

        return *this;

}

//Move assignment
lattice_site& lattice_site::operator =(lattice_site&& site1)
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
