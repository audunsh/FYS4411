#ifndef BASIS_H
#define BASIS_H

#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

class basis
{
public:
    basis(int N);
    void read(string filename, int Zn);
    void generate();
    void set_orthonormal(bool t);
    double get(int p, int q, int r, int s);
    double eval(int p, int q, int r, int s);
    double state(int p, int q, int r, int s, double D, double E);
    int Nstates;
    field<mat> v; //(Nstates, Nstates);
    field<mat> V; //(Nstates*2,Nstates*2);
    int Z;
    double h0(int i, int j);
    mat S;

};

#endif // BASIS_H
