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
    void Solve(basis Bs);

private:
<<<<<<< HEAD
    double h0(int alpha,int gamma);
    mat HF(mat C, field<mat> V);
    double calc_energy(mat C,field<mat> V);
    void get_Q();
    void get_S();
=======
    mat HF(mat C, basis Bs);
    double calc_energy(mat C,basis Bs);
>>>>>>> 3b050e3940d186edb487b3e8edede6875e365eb8
};

#endif // HFSOLVE_H
