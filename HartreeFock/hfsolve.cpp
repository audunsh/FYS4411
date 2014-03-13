#include "hfsolve.h"
#include "lib.h"
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

void HFSolve::Solve(basis Bs){
    //Solves the HF-eq's using the provided basis
    double tolerance = 10e-8;
    mat C;
    vec e_v, e_v_prev;
    Nstates = Bs.Nstates2; //set number of states in basis
    C.zeros(Nstates,Nstates);
    e_v.zeros(Nstates);
    e_v_prev.zeros(Nstates);
    for (int i = 0; i < Nstates; ++i) {
        C(i,i) = 1.0;
    }

    int iters = 0;
    for(int i=0; i<Nstates;i++){e_v_prev(i) = 1.0;} // safety margin
    while (abs(e_v.min() - e_v_prev.min()) > tolerance){ // convergence test
        iters = iters + 1;
        e_v_prev = e_v;
        // return the eigenvalues of the HF-mx to e_v and the eigenvectors to C.
        eig_sym(e_v,C,HF(C,Bs));
        C = trans(C);
        if(iters>100){
            cout << "Maximum number of iterations (100) exceeded." << endl;
            break;}
    }
    cout << "------------------------------" << endl;
    cout << "iterations: " << iters << endl;
    //cout << "eigenvalues: " << endl;
    int minIndex = 0;
    for (int i = 0; i < Nstates; ++i) {
        if (e_v[i] < e_v[minIndex]){
            minIndex = i;
        }
        //cout << e_v[i] << endl;
    }
    double E = calc_energy(C,Bs);
    cout << "Ground State Energy: " << E << endl;

}





mat HFSolve::HF(mat C, basis Bs){
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
                        interaction += C(p,beta)*C(p,delta)*Bs.V(alpha,beta)(gamma,delta);
                    }
                }
            }
            HFmx(alpha,gamma) = Bs.h0(alpha,gamma) + interaction;
        }
    }
    return HFmx;
}



double HFSolve::calc_energy(mat C, basis Bs){

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
        Energy += Bs.h0(alpha,alpha)*Ci.at(alpha,alpha);
        for(int beta = 0; beta < Nstates; beta++){
            for(int gamma = 0; gamma<Nstates; gamma++){
                for(int delta = 0; delta <Nstates; delta++){
                    Energy += 0.5 * Ci.at(alpha, gamma) * Ci.at(beta, delta) * Bs.V(alpha, beta)(gamma, delta);
                }
            }
        }
    }
    return Energy;
}




