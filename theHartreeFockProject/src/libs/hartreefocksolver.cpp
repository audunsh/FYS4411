#include "hartreefocksolver.h"
#include <basis.h>

hartreefocksolver::hartreefocksolver(){}

hartreefocksolver::hartreefocksolver(basis BS, int N, int Z){
    //initialize solver with given basis
    //This solver follows the algo described in Thijssen, p74-76

    Bs = BS;
    nStates = Bs.Nstates;        //set number of states in basis
    nElectrons = N;
    nProtons = Z;
    C.zeros(nStates,nStates);    //set initial C equal to the unit matrix
    F.zeros(nStates,nStates);    //initialize Fock matrix
    P.zeros(nStates,nStates);    //initialize Density matrix
    U.zeros(nStates,nStates);    //initialize Unitary matrix
    G.zeros(nStates,nStates);    //initialize Fock-component matrix

    Fprime.zeros(nStates,nStates); //transformed Fock matrix
    epsilon.zeros(nStates);    //eigenvalues from diagonalization
    epsilon_prev.zeros(nStates); //eigenvalues from previous diagonalization

}

double hartreefocksolver::solve(){
    //carefully following the steps laid out on pages 74-77 in Thijssen
    setupUnitMatrices();
    setupP();
    setupF();
    iterations = 0;
    while(convergenceCriteria()){
        diagonalizeF();
        updateP();
        iterations += 1;
    }
    return energy();
}

void hartreefocksolver::setupUnitMatrices(){
    //Bring overlap matrix to unit form, set up V
    eig_sym(s_diag,U,Bs.S);           //following Thijssen, p38-39
    V = U*diagmat(1.0/sqrt(s_diag));
}

void hartreefocksolver::setupP(){
    //setup density matrix, make a first guess
    P.zeros(); //we don't have any reason to do this any other ways yet.
}

void hartreefocksolver::setupF(){
    //set up the Fock matrix
    for(int p=0;p<nStates;p++){
        for(int q=0;q<nStates;q++){
            for(int r=0;r<nStates;r++){
                for(int s=0;s<nStates;s++){
                    G(p,q) += P(r,s)*(Bs.v(p,r)(q,s) - .5*Bs.v(p,r)(s,q));
                }
            }
            F(p,q) = G(p,q) + Bs.h(p,q);
        }
    }
}

void hartreefocksolver::diagonalizeF(){
    //diagonalize the Fock matrix
    Fprime = V.t()*F*V;
    eig_sym(epsilon, Cprime, Fprime);
    C = V*Cprime;
}

void hartreefocksolver::normalizeC(){
    //normalizing the C matrix, following Thijessen, p 76
    double result;
    for(int k = 0; k<nElectrons;k++){
        result = 0;
        for(int p = 0; p<nStates;p++){
            for(int q = 0; q<nStates;q++){
                result += C(p,k)*Bs.S(p,q)*C(q,k);
            }
        }
        C.col(k)/=sqrt(result);
    }
}

void hartreefocksolver::updateP(){
    //construct new density matrix
    mat Ptemp = 2*C.cols(0,nElectrons/2-1)*C.cols(0,nElectrons/2 - 1).t();
    P = 0.5*P + 0.5*Ptemp;
}

bool hartreefocksolver::convergenceCriteria(){
    //check for convergence
    bool condition = true;
    if(iterations>100){
        condition = false;
    }
    return condition;
}

double hartreefocksolver::energy(){
    //return ground state energy
    double e0 = 0;
    for(int p = 0; p<nStates;p++){
        for(int q = 0; q<nStates; q++){
            e0 += P(p,q)*Bs.h(p,q)+nProtons*Bs.nuclearPotential(p,q); //is the nuclear potential term correctly placed here?
        }
    }
    for(int p = 0; p<nStates;p++){
        for(int q = 0; q<nStates; q++){
            for(int r = 0; r<nStates;r++){
                for(int s = 0; s<nStates; s++){
                    e0 += 0.25*P(p,q)*P(s,r)*(2*Bs.v(p,q)(r,s)-Bs.v(p,q)(s,r)); //need to ask about this
                }
            }
        }
    }
    return e0;
}


