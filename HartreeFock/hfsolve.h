#ifndef HFSOLVE_H
#define HFSOLVE_H
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

class HFSolve{
    int Z,N;
public:
    HFSolve(int Zn, int Nn);

    field<mat> init(string filename);
    double state(int p, int q, int r, int s, double D, double Ex);
    double return_init();
    void Solve(field<mat> V);

private:
    double h0(int alpha,int gamma);
    mat HF(mat C, field<mat> V);
    double calc_energy(mat C,field<mat> V);
};

#endif // HFSOLVE_H