#ifndef HFSOLVE_H
#define HFSOLVE_H
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

class HFSolve{
    int Z,N, Nstates;
public:
    HFSolve(int Zn, int Nn, int Ns);

    field<mat> init(string filename);

    void use_basis(string choise_of_basis); // GTO,STO,filename

    double state(int p, int q, int r, int s, double D, double Ex);

    double return_init();

    void Solve(field<mat> V);

private:
    double h0(int alpha,int gamma);
    mat HF(mat C, field<mat> V);
    double calc_energy(mat C,field<mat> V);
    void get_Q();
    void get_S();
};

#endif // HFSOLVE_H
