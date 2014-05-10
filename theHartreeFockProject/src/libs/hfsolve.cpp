#include <hfsolve.h>
#include <lib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <basis.h>


using namespace arma;

HFSolve::HFSolve(int Zn, int Nn){
    Z = Zn;
    N = Nn;
}

void HFSolve::init_solver(){
    //initialize and set up solver when basis is supplied
    Nstates = Bs.Nstates;        //set number of states in basis
    C.zeros(Nstates,Nstates);    //set initial C equal to the unit matrix
    F.zeros(Nstates,Nstates);    //initialize Fock matrix
    P.zeros(Nstates,Nstates);    //initialize Density matrix
    for (int i = 0; i < Nstates; ++i) {
        C(i,i) = 1.0;
    }
    e_v.zeros(Nstates);          //vec containing current eigenvalues of transformed Fock-matrix
    e_v_prev.zeros(Nstates);     //vec containing previous eigenvalues of transformed Fock-matrix
    eig_sym(s,U,Bs.S);           //set up the transformed equation following Thijssen, p38-39
    V = U*diagmat(1.0/sqrt(s));
    C_trans = inv(V)*C;
    setupP(C);                   //set up the density matrix

}

void HFSolve::updateF(){
    //Update the Fock matrix using P and basis Bs
    for(int p = 0; p<Nstates;p++){
        for(int q = 0; q<Nstates;q++){
            F(p,q) = Bs.h(p,q);
            for(int r = 0; r<Nstates;r++){
                for(int s = 0; s<Nstates;s++){
                    F(p,q) += 0.5 * P(r,s)* (2*Bs.v(p,q)(r,s) - Bs.v(p,q)(s,r));
                }
            }
        }
    }
}

void HFSolve::advance(){
    //For each iteration in the HF-algorithm
    e_v_prev = e_v;
    updateF();
    F_trans = V.t()*F*V; //transforming the hartree fock matrix
    // return the eigenvalues of the HF-mx to e_v and the eigenvectors to C.
    eig_sym(e_v,C_trans,F_trans); //solving the transformed equation
    C = V*C_trans;    //transforming back
    C = trans(C);     //transposing C
    normalize_col(C); //normalizing the coloumns of C
    setupP(C);

}

double HFSolve::Solve(basis BS){
    //Solves the HF-eq's using the provided basis
    double tolerance = 10e-8;


    Bs = BS;
    init_solver();

    int iters = 0;
    for(int i=0; i<Nstates;i++){e_v_prev(i) = 1.0;}      //Set convergence parameter to false, see while-loop below
    while (abs(e_v.min() - e_v_prev.min()) > tolerance){ // convergence test
        iters = iters + 1;

        advance();


        if(iters>100){
            cout << "Maximum number of iterations (100) exceeded." << endl;
            break;}
    }


    cout << "------------------------------" << endl;
    cout << "iterations: " << iters << endl;
    int minIndex = 0;
    for (int i = 0; i < Nstates; ++i) {
        if (e_v[i] < e_v[minIndex]){
            minIndex = i;
        }
        //cout << e_v[i] << endl;
    }
    C.print();
    double E = calc_energy(C);
    //E = energy();

    //cout << "Ground State Energy: " << E << endl;
    return E;
}





mat HFSolve::HFmatrix(mat C){
    /* Sets up the Hartree-Fock matrix
     * using the coefficients given in C
     */

    mat HFmx;
    HFmx.zeros(C.n_rows, C.n_cols);
    for (int alpha = 0; alpha < Nstates; ++alpha) {
        for (int gamma = 0; gamma < Nstates; ++gamma) {
            double interaction = 0;
            for (int p = 0; p < N; ++p) {
                for (int beta = 0; beta < Nstates; ++beta) {
                    for (int delta = 0; delta < Nstates; ++delta) {
                        interaction += 0.5*C(p,beta)*C(p,delta)*(2*Bs.v(alpha,beta)(gamma,delta)-Bs.v(alpha,beta)(delta,gamma));
                    }
                }
            }
            HFmx(alpha,gamma) = Bs.h(alpha,gamma) + interaction;
        }
    }
    return HFmx;
}



double HFSolve::calc_energy(mat C){

    mat Ci = zeros(Nstates,Nstates);
    double s;
    for(int p= 0; p<Nstates;p++){
        for(int q=0; q<Nstates; q++){
            s = 0;
            for(int k=0; k<N; k++){
                s+=C.at(k,p) * C.at(k,q);
                Ci.at(p,q) = s;
            }
        }
    }
    double Energy = 0;
    for(int alpha = 0; alpha < Nstates; alpha++){
        Energy += Bs.h(alpha,alpha)*Ci.at(alpha,alpha);
        for(int beta = 0; beta < Nstates; beta++){
            for(int gamma = 0; gamma<Nstates; gamma++){
                for(int delta = 0; delta <Nstates; delta++){
                    //Energy += 0.5 * Ci.at(alpha, gamma) * Ci.at(beta, delta) * Bs.v(alpha, beta)(gamma, delta);
                    Energy += 0.5 * Ci.at(alpha, gamma) * Ci.at(beta, delta) * (2*Bs.v(alpha, beta)(gamma, delta)-Bs.v(alpha, beta)(delta,gamma));
                }
            }
        }
    }
    return Energy;
}

void HFSolve::normalize_col(mat C){
    //normalizing the C matrix, following Thijessen, p 76
    double result;
    for(int k = 0; k<Nstates;k++){
        result = 0;
        for(int p = 0; p<Nstates;p++){
            for(int q = 0; q<Nstates;q++){
                result += C(p,k)*Bs.S(p,q)*C(q,k);
            }
        }
        C.col(k)/=sqrt(result);
    }
}

void HFSolve::setupP(mat C){
    mat P = 2*C*C.t();
}

double HFSolve::energy(){
    return 0.5*accu(P % (Bs.h + F));//+Bs.nuclearPotential;
}
