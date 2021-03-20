#include <iostream>
#include <complex>
#include <cmath>
#include <random>
#include <chrono>
#include "rand.h"
//#include "lattice.h"

unsigned seedGen()
{
        typedef std::chrono::high_resolution_clock myclock;
        myclock::time_point beginning = myclock::now();
        myclock::duration d = myclock::now() - beginning;
        unsigned seed2 = d.count();
        return seed2;
}


const matrix_complex normalHermitianMatrix( std::mt19937_64& g ) //Also traceless
{
        double a,b,c;
        std::normal_distribution<double> distribution(0.0,1.0);
        matrix_complex herm;
        a = distribution(g);
        b = distribution(g);
        c = distribution(g);
        herm(0,0) = complex<double>(c);
        herm(1,1) = complex<double>(-c);
        herm(0,1) = complex<double>(a,-b);
        herm(1,0) = complex<double>(a,b);
        return herm;
}

const matrix_complex uniformSU2Matrix( std::mt19937_64& g ) // Identity plus an anti-Hermitian matrix
{
        double a[4], total;
        matrix_complex su;
        total = 0;
        for(int i=0; i<4; i++)
        {
                a[i] = uniformReal(g,-1.0,1.0);
                total += a[i] * a[i];
        }
        total = 1/sqrt(total);
        for(int i=0; i<4; i++)
                a[i] *= total;
        su(0,0) = complex<double>(a[0],a[1]);
        su(1,1) = complex<double>(a[0],-a[1]);
        su(0,1) = complex<double>(a[2],a[3]);
        su(1,0) = complex<double>(-a[2],a[3]);
        return su;
}

const matrix_complex smallSU2Matrix( std::mt19937_64& g )
{
        double a[4], total,epsilon;
        matrix_complex su;
        total = 0;
        for(int i=1; i<4; i++)
        {
                a[i] = uniformReal(g,-1.0,1.0);
                total += a[i] * a[i];
        }
        epsilon = uniformReal(g,0.0,EPSILON_MAX);
        total = epsilon/sqrt(total);
        for(int i=1; i<4; i++)
                a[i] *= total;
        a[0] = sqrt(1.0d - epsilon*epsilon);
        su(0,0) = complex<double>(a[0],a[1]);
        su(1,1) = complex<double>(a[0],-a[1]);
        su(0,1) = complex<double>(a[2],a[3]);
        su(1,0) = complex<double>(-a[2],a[3]);
        return su;
}

const matrix_complex smallHermitianMatrix( std::mt19937_64& g )
{
        double a,b,c;

        matrix_complex herm;
        a = HERM_STEP_MAX * uniformReal(g,-1.0,1.0);
        b = HERM_STEP_MAX * uniformReal(g,-1.0,1.0);
        c = HERM_STEP_MAX * uniformReal(g,-1.0,1.0);
        herm(0,0) = complex<double>(c);
        herm(1,1) = complex<double>(-c);
        herm(0,1) = complex<double>(a,-b);
        herm(1,0) = complex<double>(a,b);
        return herm;
}

const double uniformReal(std::mt19937_64& g)
{
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        return distribution(g);
}

const double uniformReal(std::mt19937_64& g, double low, double high)
{
        std::uniform_real_distribution<double> distribution(low,high);
        return distribution(g);
}
