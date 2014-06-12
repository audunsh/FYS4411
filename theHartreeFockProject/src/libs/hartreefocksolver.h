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
    void setupP();       //setup density matrix, make a first guess
    void setupF();       //setup the Fock matrix
    void diagonalizeF(); //diagonalize the Fock matrix
    void normalizeC();   //normalize the coefficient matrix
    void updateP();      //construct new density matrix
    void setupCoupledMatrix(); //set up direct and exchange terms from basis
    void printMatrices();//print the iterating matrices
    bool convergenceCriteria(); //check for convergence
    double energy();     //return ground state energy
    double solve();      //automated solving process, returns energy
    double coupledMatrixTilde(int p, int q, int r, int s); //return direct and exchange term

    //new code to debug, possibly for deletion
    void setupCoupledMatrix_unused();
    double calcEnergy2();
    double energyCalc();
    double evaluateProbabilityDensity(vec3 r);
    void createDensityMap(string filename);
    void reset(basis BS, int N, int Z);

private:
    basis Bs;

    field<mat> coupledMatrix;

    cube densityMap;

    mat V;      //Transformation matrix
    mat P;      //Density matrix
    mat P_prev; //previous density matrix
    mat F;      //Fock matrix
    mat C;      //Coefficient matric
    mat Fprime; //transformed Fock matrix
    mat Cprime; //transformed Coefficient matric

    mat G;      //Coulomb and exchange contribution matrix
    mat U;      //unitary matrix

    vec epsilon;      //eigenvalues from current diagonalization
    vec epsilon_prev; //eigenvalues from previous diagonalization
    vec s_diag;       //diagonalized S matrix

    int nElectrons; //number of electrons
    int nStates;    //number of states
    int nProtons;   //number of protons, can be removed
    int iterations; //number of iterations (counter)

    double energyPrev = 10e10;
    double tolerance = 10e-10;
    double dampingFactor = 0.95;
};

#endif // HARTREEFOCKSOLVER_H
