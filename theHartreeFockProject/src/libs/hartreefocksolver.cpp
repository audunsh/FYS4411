#include "hartreefocksolver.h"
#include <basis.h>

hartreefocksolver::hartreefocksolver(){}

hartreefocksolver::hartreefocksolver(basis BS, int N, int Z){
    //initialize solver with given basis
    //This solver follows the algo described in Thijssen, p74-76
    Bs = BS;
    nStates = Bs.Nstates;        //set number of states in basis
    nElectrons = N;              //set number of electrons
    nProtons = Z;                //set number of protons
    //initializing all matrices and vectors
    C.zeros(nStates,nElectrons/2);    //set initial C equal to the unit matrix
    F.zeros(nStates,nStates);    //initialize Fock matrix
    P.zeros(nStates,nStates);    //initialize Density matrix
    U.zeros(nStates,nStates);    //initialize Unitary matrix
    G.zeros(nStates,nStates);    //initialize Fock-component matrix
    Fprime.zeros(nStates,nStates);//transformed Fock matrix
    epsilon.zeros(nStates);       //eigenvalues from diagonalization
    epsilon_prev.zeros(nStates);  //eigenvalues from previous diagonalization
    setupCoupledMatrix();
    setupP();
    setupF();
    s_diag.zeros(nStates);
}

double hartreefocksolver::solve(){
    //carefully following the steps laid out on pages 74-77 in Thijssen
    setupUnitMatrices();
    setupP();
    iterations = 0;
    //printMatrices();
    while(convergenceCriteria()){
        epsilon_prev = epsilon;
        energyPrev = energy();

        setupF();
        diagonalizeF();
        normalizeC();
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
    //Comment, 10/5/14
    /* I'm not sure about the steps laid out by Thijssen at page 75; "Calculate Coulomb and exchange contribution ..."
     * How does the matrix G below compare to Thijssens? (G as the sum over r,s)
     * Update: Changed indices to match www.github.com/dragly/hartree-fock
     */
    for(int p=0;p<nStates;p++){
        for(int q=0;q<nStates;q++){
            F(p,q) = Bs.h(p,q);
            for(int r=0;r<nStates;r++){
                for(int s=0;s<nStates;s++){
                    //F(p,q) += 0.5*P(s,r)*coupledMatrixTilde(p,q,r,s);
                    F(p,q) += P(s,r) * (coupledMatrix(p,r)(q,s)-0.5*coupledMatrix(p,r)(s,q));
                }
            }
        }
    }
}

void hartreefocksolver::diagonalizeF(){
    //diagonalize the Fock matrix
    Fprime = V.t()*F*V;
    eig_sym(epsilon, Cprime, Fprime);
    C = V*Cprime.submat(0, 0, nStates - 1, nElectrons/2 -1);
}

void hartreefocksolver::normalizeC(){
    //normalizing the C matrix, following Thijessen, p 76
    double result;
    for(int k = 0; k<nElectrons/2;k++){
        result = 0.0;
        for(int p = 0; p<nStates;p++){
            for(int q = 0; q<nStates;q++){
                result += C(p,k)*Bs.S(p,q)*C(q,k);
            }
        }
        C.col(k)=C.col(k)/sqrt(result);
    }
}

void hartreefocksolver::updateP(){
    //construct new density matrix
    //following Dragly at github.com/dragly/hartree-fock
    mat Ptemp = 2*C*C.t();
    if(P.n_elem>0){
        P = .5 * P + (1 - .5)* Ptemp;
    }
    else{
        P = Ptemp;
    }
}

bool hartreefocksolver::convergenceCriteria(){
    //Evaluate convergence conditions
    bool condition = true;
    if(iterations>100){
        condition = false;
    }
    if(abs(energyPrev-energy())<tolerance){
        condition = false;
    }
    return condition;
}

double hartreefocksolver::energy(){
    //return ground state energy, following Svenn-Arne Dragly at
    //www.github.com/dragly/hartree-fock
    double e0 = 0;
    for(int p = 0; p<nStates;p++){
        for(int q = 0; q<nStates; q++){
            e0 += P(p,q)*Bs.h(p,q);
            for(int r = 0; r<nStates;r++){
                for(int s = 0; s<nStates; s++){
                    e0 += 0.25*P(p,q)*P(s,r)*coupledMatrixTilde(p,q,r,s);
                    //e0 += 0.5*P(p,q)*P(s,r)*(coupledMatrix(p,r)(q,s)-0.5*coupledMatrix(p,r)(s,q));

                }
            }
        }
    }
    return e0+Bs.nnInteraction();
}

double hartreefocksolver::coupledMatrixTilde(int p, int q, int r, int s){
    //return direct and exchange term, following Svenn-Arne Dragly at
    //www.github.com/dragly/hartree-fock
    return 2*coupledMatrix(p,r)(q,s) - coupledMatrix(p,r)(s,q);

}

void hartreefocksolver::setupCoupledMatrix_unused(){
    //still following Dragly, further references to indexation differences between Thijssen and Helgaker
    coupledMatrix.set_size(nStates,nStates);
    for(int p = 0; p<nStates;p++){
        for(int r = 0; r<nStates;r++){
            coupledMatrix(p,r)=zeros(nStates,nStates);
            for(int q = p; q<nStates;q++){
                for(int s = r; s<nStates;s++){
                    coupledMatrix(p,r)(q,s) = Bs.v(p,q)(r,s);
                }
            }
        }
    }
    double val;
    for(int p = 0; p<nStates;p++){
        for(int r = 0; r<nStates;r++){
            for(int q = p; q<nStates;q++){
                for(int s = r; s<nStates;s++){
                    val = coupledMatrix(p,r)(q,s);
                    coupledMatrix(q,s)(p,r) = val;
                    coupledMatrix(q,r)(p,s) = val;
                    coupledMatrix(p,s)(q,r) = val;
                    coupledMatrix(r,p)(s,q) = val;
                    coupledMatrix(s,p)(r,q) = val;
                    coupledMatrix(r,q)(s,p) = val;
                    coupledMatrix(s,q)(r,p) = val;
                }
            }
        }
    }
    coupledMatrix.print();
}

void hartreefocksolver::printMatrices(){
    cout << "Fock matrix" << endl;
    F.print();
    cout << " " << endl;

    cout << "Coeff matrix" << endl;
    C.print();
    cout << " " << endl;

    cout << "Density matrix" << endl;
    P.print();
    cout << " " << endl;
    cout << "----------------------" << endl;

}

void hartreefocksolver::setupCoupledMatrix(){
    int n = Bs.Nstates;
    coupledMatrix.set_size(n, n);
    for (int p = 0; p<n; p++){
        for (int q = 0; q<n; q++){
            coupledMatrix(p, q) = zeros(n, n);
        }
    }
    for (int p = 0; p<n; p++){
        for (int q = 0; q<n; q++){
            for (int r = 0; r<n; r++){
                for (int s = 0; s<n; s++){
                    coupledMatrix(p, r)(q, s) = Bs.v(p, r)(q, s);
                }
            }
        }
    }
}
