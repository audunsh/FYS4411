#include "hfsolve.h"
#include "lib.h"
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <basis.h>

using namespace arma;

HFSolve::HFSolve(int Zn, int Nn, int Ns){
    Z = Zn;
    N = Nn;
    Nstates = Ns;
}
field<mat> HFSolve::init(string filename){
    /* Read filename, and set up initial matrix
     * Z is the atomic number
     * N is the number of electrons
     */
    //To be implemented: This part of the solver can be migrated to a separate "basis"-class.
    ifstream myfile;
    myfile.open(filename.c_str());
    field<mat> v(3,3);     // 3x3 field matrix with matrix elements
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            v(i,j) = zeros<mat>(3,3); // fill V with 3x3 mx elements
        }
    }

    if (myfile.is_open()){
        int p,q,r,s;
        double value;
        while (!myfile.eof()){
            myfile >> p;
            myfile >> q;
            myfile >> r;
            myfile >> s;
            myfile >> value;
            v(p,q)(r,s) = value;
        }
    }
    else
        cout << "Did not manage to open file in HFSolve::init()"<< endl;

    field<mat> V(6,6);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            V(i,j) = zeros<mat>(6,6); // fill V with 3x3 mx elements
        }
    }
    double D = 0;
    double Ex = 0;
    for (int p = 0; p < Nstates; ++p) {
        for (int q = 0; q < Nstates; ++q) {
            for (int r = 0; r < Nstates; ++r) {
                for (int s= 0; s < Nstates; ++s) {
                    D = v(p/2,q/2)(r/2,s/2);  // Direct term
                    Ex = v(p/2,q/2)(s/2,r/2); // Exchange term
                    V(p,q)(r,s) = Z*state(p,q,r,s,D,Ex);
                }
            }
        }
    }

    return V;
}

double HFSolve::state(int p, int q, int r, int s, double D, double Ex){
    //Evaluating spin configuration, returning direct and/or exchange term or 0
    double S;
    int s1,s2,s3,s4;
    s1 = p%2;
    s2 = q%2;
    s3 = r%2;
    s4 = s%2;
    if (s1 == s2){
        if (s3 == s4){
            if ( s1 == s3){
                S = D-Ex;
            }
            else{
                S = 0;
            }
        }
    }
    else if (s1 != s2){
        if (s3 != s4){
            if (s1 == s3){
                S = D;
            }
            else{
                S = -Ex;
            }
        }
        else{
            S = 0;
        }
    }
    return S;
}

void HFSolve::get_Q(){
    /* Get the basis
     */
}

void HFSolve::get_S(){

}

double HFSolve::h0(int alpha,int gamma){
    // the one-body interaction
    double h = 0;
    if (alpha == gamma){
        double n = alpha/2 + 1.0;
        h = -(Z*Z)/(2*n*n);
    }
    return h;
}



mat HFSolve::HF(mat C, field<mat> V){
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
                        interaction += C(p,beta)*C(p,delta)*V(alpha,beta)(gamma,delta);
                    }
                }
            }
            HFmx(alpha,gamma) = h0(alpha,gamma) + interaction;
        }
    }
    return HFmx;
}

void HFSolve::Solve(field<mat> V){
    /* Sets up the identity matrix as the initial guess
     * on how the coefficient matrix C should look like
     * And solves the HF equations
     */
    double tolerance = 10e-8;
    mat C;
    vec e_v, e_v_prev;
    C.zeros(6,6);
    e_v.zeros(6);
    e_v_prev.zeros(6);
    for (int i = 0; i < 6; ++i) {
        C(i,i) = 1.0;
    }

    int iters = 0;
    for(int i=0; i<6;i++){e_v_prev(i) = 1.0;} // safety margin
    while (abs(e_v.min() - e_v_prev.min()) > tolerance){ // convergence test
        iters = iters + 1;
        e_v_prev = e_v;
        // return the eigenvalues of the HF-mx to e_v and the eigenvectors to C.
        eig_sym(e_v,C,HF(C,V));
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
    double E = calc_energy(C,V);
    cout << "Ground State Energy: " << E << endl;


}

double HFSolve::calc_energy(mat C, field<mat> V){

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
        Energy += h0(alpha,alpha)*Ci.at(alpha,alpha);
        for(int beta = 0; beta < Nstates; beta++){
            for(int gamma = 0; gamma<Nstates; gamma++){
                for(int delta = 0; delta <Nstates; delta++){
                    Energy += 0.5 * Ci.at(alpha, gamma) * Ci.at(beta, delta) * V(alpha, beta)(gamma, delta);
                }
            }
        }
    }
    return Energy;
}




