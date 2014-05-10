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
    double energy();
    void init_solver();
    void advance();


private:
    void updateF();
    void normalize_col(mat C);
    void setupP(mat C);
    void setupG();



    mat HFmatrix(mat C);
    double calc_energy(mat C);

    mat G,P,C,U, V, F,F_trans, C_trans; //The transformed matrices (Thijessen, p38-39)
    vec s, e_v, e_v_prev;
    basis Bs;

};

#endif // HFSOLVE_H
