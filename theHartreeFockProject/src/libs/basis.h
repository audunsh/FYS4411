#ifndef BASIS_H
#define BASIS_H

#include <string>
#include <armadillo>

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
    double get(int p, int q, int r, int s);                       //function to retrieve two-body integral (precalculated or otherwise
    double eval(int p, int q, int r, int s);                      //function to evaluate two-body integral
    double state(int p, int q, int r, int s, double D, double E); //function to evaluate spin-dependence
    int Nstates, Nstates2;                                        //number of basis states (spin not included)
    field<mat> v;                                                 //basis, spin not included
    field<mat> V;                                                 //basis with spin included
    int Z;                                                        //interaction parameter (taken to be number of protons)
    double h0(int i, int j);                                      //one body interaction
    mat S,h,H;                                                     //overlap matrix, onebody interaction matrices
};

#endif // BASIS_H
