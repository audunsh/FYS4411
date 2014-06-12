#include "hartreefocksolver.h"
#include <basis.h>

hartreefocksolver::hartreefocksolver(){}

hartreefocksolver::hartreefocksolver(basis BS, int N, int Z){
    //initialize solver with given basis, number of electrons (and now superfluos number of protons)
    //This solver follows the algo described in Thijssen, p74-76
    Bs = BS;
    nStates = Bs.Nstates;        //set number of states in basis
    nElectrons = N;              //set number of electrons
    nProtons = Z;                //set number of protons, may be removed

    //initializing all matrices and vectors
    C.zeros(nStates,nElectrons/2);//set initial C equal to the unit matrix
    F.zeros(nStates,nStates);     //initialize Fock matrix
    P.zeros(nStates,nStates);     //initialize Density matrix
    U.zeros(nStates,nStates);     //initialize Unitary matrix
    G.zeros(nStates,nStates);     //initialize Fock-component matrix
    Fprime.zeros(nStates,nStates);//transformed Fock matrix
    epsilon.zeros(nStates);       //eigenvalues from diagonalization
    epsilon_prev.zeros(nStates);  //eigenvalues from previous diagonalization

    setupCoupledMatrix();  //import particle-particle interaction integrals
    //setupCoupledMatrix_unused();  //import particle-particle interaction integrals
    setupP();              //set up initial density matrix
    s_diag.zeros(nStates);
}

void hartreefocksolver::reset(basis BS, int N, int Z){
    //initialize solver with given basis, number of electrons (and now superfluos number of protons)
    //This solver follows the algo described in Thijssen, p74-76
    Bs = BS;
    nStates = Bs.Nstates;        //set number of states in basis
    nElectrons = N;              //set number of electrons
    nProtons = Z;                //set number of protons, may be removed

    //initializing all matrices and vectors
    C.zeros(nStates,nElectrons/2);//set initial C equal to the unit matrix
    F.zeros(nStates,nStates);     //initialize Fock matrix
    P.zeros(nStates,nStates);     //initialize Density matrix
    U.zeros(nStates,nStates);     //initialize Unitary matrix
    G.zeros(nStates,nStates);     //initialize Fock-component matrix
    Fprime.zeros(nStates,nStates);//transformed Fock matrix
    epsilon.zeros(nStates);       //eigenvalues from diagonalization
    epsilon_prev.zeros(nStates);  //eigenvalues from previous diagonalization

    setupCoupledMatrix();  //import particle-particle interaction integrals
    setupP();              //set up initial density matrix
    s_diag.zeros(nStates);
}

double hartreefocksolver::solve(){
    //carefully following the steps laid out on pages 74-77 in Thijssen

    setupUnitMatrices();
    setupP();
    iterations = 0;

    setupF();
    //printMatrices();
    while(convergenceCriteria()){
        epsilon_prev = epsilon;
        energyPrev = energyCalc();

        setupF();

        diagonalizeF();
        normalizeC();
        updateP();

        iterations += 1;
        //cout << energy() << endl;
    }
    //cout << "Converged in " << iterations << " iterations." << endl;
    //printMatrices();
    //createDensityMap();
    //return iterations;
    return energyCalc();
}

double hartreefocksolver::evaluateProbabilityDensity(vec3 r){
    double result = 0;
    for(int p=0;p<nStates;p++){
        for(int q=0;q<nStates;q++){
            result += P(p,q)*Bs.evaluateContracted(p, r)*Bs.evaluateContracted(q,r);
        }
    }
    return result;
}

void hartreefocksolver::createDensityMap(string filename){
    int dim=1000;
    double dx = 0.01;
    double dV = dx*dx*dx;
    densityMap.zeros(dim,dim,dim);
    mat densitySlice;
    densitySlice.zeros(dim,dim);
    vec3 R;
    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            R = {(double) i*dx,(double) j*dx, 0.0};
            densitySlice(i,j)=evaluateProbabilityDensity(R);
            /*
            for(int k=-dim;k<dim;k++){
                R = {(double) i*dx,(double) j*dx, (double) k*dx};
                densityMap(i+dim,j+dim,k+dim) = evaluateProbabilityDensity(R)*dV;
            */
        }
    }
    //densityMap.print();
    //densityMap.save("testmap3", raw_ascii);
    //densitySlice.print();
    densitySlice.save(filename, raw_ascii);

}

void hartreefocksolver::setupCoupledMatrix(){
    //import particle-particle interaction integrals
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
                    //(p, r)(q, s) = Bs.v(p, r)(q, s); //alt (1) "Strange" (ALT 1), working
                    coupledMatrix(p, r)(q, s) = Bs.v(p, q)(r, s); //alt (1) "Strange" (CHANGED TODAY)

                }
            }
        }
    }
}

void hartreefocksolver::setupF(){
    //set up the Fock matrix
    for(int p=0;p<nStates;p++){
        for(int q=0;q<nStates;q++){
            F(p,q) = Bs.h(p,q);
            for(int r=0;r<nStates;r++){
                for(int s=0;s<nStates;s++){
                    //F(p,q) += 0.5*coupledMatrixTilde(p,q,r,s)*P(s,r);  //Alt. 1 "Thijssen"
                    F(p,q) += 0.5*coupledMatrixTilde(p,q,r,s)*P(r,s);  //Alt. 2 27/5 2014
                }
            }
        }
    }
}

double hartreefocksolver::energyCalc(){
    //return energy for current Fock matrix
    return 0.5*accu(P % (Bs.h + F))+Bs.nnInteraction();
}

double hartreefocksolver::energy(){
    //return ground state energy
    double e0 = 0;
    for(int p = 0; p<nStates;p++){
        for(int q = 0; q<nStates; q++){
            e0 += P(p,q)*Bs.h(p,q);
            for(int r = 0; r<nStates;r++){
                for(int s = 0; s<nStates; s++){
                    e0 += 0.25*coupledMatrixTilde(p,q,r,s)*P(p,q)*P(s,r);
                }
            }
        }
    }
    return e0+Bs.nnInteraction();
}

double hartreefocksolver::coupledMatrixTilde(int p, int q, int r, int s){
    //return direct and exchange term, weigthed to include spin
    //return 2*coupledMatrix(p,q)(r,s) - coupledMatrix(p,s)(r,q);
    //return 2*coupledMatrix(p,r)(q,s) - coupledMatrix(p,r)(s,q); //WORKING, FAMILY A
    return 2*coupledMatrix(p,r)(q,s) - coupledMatrix(p,r)(s,q); //CHANGED 27/5
}

void hartreefocksolver::setupUnitMatrices(){
    //Bring overlap matrix to unit form, set up V
    eig_sym(s_diag,U,Bs.S);           //following Thijssen, p38-39
    V = U*diagmat(1.0/sqrt(s_diag));
}

void hartreefocksolver::setupP(){
    //setup density matrix, make a first guess
    //P.zeros(); //we don't have any reason to do this any other ways yet.
    P.eye();
}

void hartreefocksolver::diagonalizeF(){
    //diagonalize the Fock matrix
    Fprime = V.t()*F*V;
    eig_sym(epsilon, Cprime, Fprime);
    C = V*Cprime.submat(0, 0, nStates - 1, nElectrons/2 -1);
    //C = V*Cprime;
}

void hartreefocksolver::normalizeC(){
    //Normalize the coefficient matrix
    double norm;
    for(int i = 0; i<nElectrons/2;i++){
        norm = dot(C.col(i),Bs.S*C.col(i));
        C.col(i) = C.col(i)/sqrt(norm);
    }
}

void hartreefocksolver::updateP(){
    //construct new density matrix
    //P = dampingFactor*P + (1-dampingFactor)*2.0*C.cols(0, nElectrons/2.0 - 1)*C.cols(0, nElectrons/2.0 - 1).t();
    //P =C.cols(0, nElectrons/2.0 - 1)*C.cols(0, nElectrons/2.0 - 1).t();
    //P = C*C.t();
    P = dampingFactor*P + (1-dampingFactor)*2.0*C.cols(0, nElectrons/2.0 -1)*C.cols(0, nElectrons/2.0 -1).t();
}

bool hartreefocksolver::convergenceCriteria(){
    //Evaluate convergence conditions
    bool condition = true;
    if(iterations>5000){
        condition = false;
    }
    if(abs(energyPrev-energy())<tolerance){
        condition = false;
    }
    return condition;
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

    cout << "H matrix" << endl;
    Bs.h.print();
    cout << " " << endl;
    cout << "----------------------" << endl;


    cout << "Coupled matrix" << endl;
    coupledMatrix.print();
    cout << " " << endl;
    cout << "----------------------" << endl;

    cout << "Overlap" << endl;
    Bs.S.print();
    cout << " " << endl;
    cout << "----------------------" << endl;



}


