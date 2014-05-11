#ifndef BASIS_H
#define BASIS_H

#include <string>
#include <armadillo>
#include <contracted.h>
#include <primitive.h>
#include <integrator.h>
#include <boysfunction.h>


using namespace std;
using namespace arma;

class basis
{
public:
    basis();
    void set_size(int N);                                         //set size of basis, N=number of orbitals
    void read(string filename, int Zn);                           //read basis from file
    void set_orthonormal();                                       //if true: set the overlap matrix to I
    void expand();                                                //expand basis for explicit spin-dependence
    void init_integrals();                                        //set up and solve all integrals for a gaussian basis
    void init_STO_3G(string configuration,double nProtons);                       //initialize STO-3G basis set for given configuration
    void init_HTO4(int nProtons);                                //initialize hydrogen type basis with 4 orbitals
    double state(int p, int q, int r, int s, double D, double E); //function to evaluate spin-dependence
    int Nstates, Nstates2;                                        //number of basis states (spin not included)
    field<mat> v;                                                 //basis, spin not included
    field<mat> V;                                                 //basis with spin included
    double Z;                                                     //interaction parameter (taken to be number of protons)
    double h0(int i, int j);                                      //one body interaction
    mat S,h,H,nuclearPotential;                                   //overlap matrix, onebody interaction matrices


private:
    contracted basisSet[]; //use one of these...
    vector<contracted> basisSts; //use one of these...
    int Nprimitives;

};

#endif // BASIS_H
