#ifndef HARTREEFOCKSOLVER_H
#define HARTREEFOCKSOLVER_H

#include <basis.h>
#include <armadillo>

class hartreefocksolver
{
public:
    hartreefocksolver();
    hartreefocksolver(basis BS, int N, int Z);
    void setupUnitMatrices(); //Bring overlap matrix to unit form
    void setupP(); //setup density matrix, make a first guess
    void setupF(); //setup the Fock matrix
    void diagonalizeF(); //diagonalize the Fock matrix
    void normalizeC();
    void updateP(); //construct new density matrix
    bool convergenceCriteria(); //check for convergence
    double energy(); //return ground state energy
    double solve(); //automated solving process, returns energy

private:
    basis Bs;
    mat V; //Transformation matrix
    mat P; //Density matrix
    mat P_prev; //previous density matrix
    mat F; //Fock matrix
    mat C; //Coefficient matric
    mat Fprime; //transformed Fock matrix
    mat Cprime; //transformed Coefficient matric

    mat G; //Coulomb and exchange contribution matrix
    mat U; //unitary matrix

    vec epsilon;
    vec epsilon_prev;
    vec s_diag; //diagonalized S matrix

    int nElectrons;
    int nStates;
    int nProtons;
    int iterations;

};

#endif // HARTREEFOCKSOLVER_H
