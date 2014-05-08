#ifndef HFSOLVE_H
#define HFSOLVE_H
#include <armadillo>
#include <string>
#include <basis.h>


using namespace std;
using namespace arma;

class HFSolve{
    int Z,N, Nstates;
public:
    HFSolve(int Zn, int Nn);
    double Solve(basis BS);

private:
    mat HF(mat C);
    double calc_energy(mat C);
    mat C,U, V, F_trans, C_trans; //The transformed matrices (Thijessen, p38-39)
    vec s;
    void normalize_C(mat C);
    basis Bs;

};

#endif // HFSOLVE_H
