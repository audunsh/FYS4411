#ifndef BASIS_H
#define BASIS_H

#include <string>
#include <armadillo>
#include <contracted.h>
#include <primitive.h>
#include <integrator.h>

using namespace std;
using namespace arma;

class basis
{
public:
    basis(int N);                                                 //initialize object
    void read(string filename, int Zn);                           //read basis from file
    void generate();                                              //generate basis
    void set_orthonormal();                                       //if true: set the overlap matrix to I
    void expand(); //expand basis for explicit spin-dependence
    void init_overlap();
    void init_integrals();
    void init_STO_3G(string configuration); //set up STO-3G basis set for given configuration
    double get(int p, int q, int r, int s);                       //function to retrieve two-body integral (precalculated or otherwise
    double eval(int p, int q, int r, int s);                      //function to evaluate two-body integral
    double state(int p, int q, int r, int s, double D, double E); //function to evaluate spin-dependence
    int Nstates, Nstates2;                                        //number of basis states (spin not included)
    field<mat> v;                                                 //basis, spin not included
    field<mat> V;                                                 //basis with spin included
    int Z;                                                        //interaction parameter (taken to be number of protons)
    double h0(int i, int j);                                      //one body interaction
    mat S,h,H;                                                     //overlap matrix, onebody interaction matrices
    //we also need an array of contracted containing the primitives in each orbital
    //this will constitute the basis
private:
    contracted basisSet[];
    int Nprimitives;

};

#endif // BASIS_H
