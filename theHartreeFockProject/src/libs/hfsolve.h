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
    double Solve(basis Bs);

private:
    mat HF(mat C, basis Bs);
    double calc_energy(mat C,basis Bs);
};

#endif // HFSOLVE_H
