#include "hfsolve.h"
#include "lib.h"
#include <fstream>
#include <string>
#include <iostream>
#include <vector>

using namespace arma;

HFSolve::HFSolve(int Zn, int Nn){
    Z = Zn;
    N=Nn;
}
field<mat> HFSolve::init(string filename){
    /* Read filename, and set up initial matrix
     * Z is the atomic number
     * N is the number of electrons
     */
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
            //cout << p << " " << q << " " << r << " " << s << " " << value << " " << v(p,q)(r,s) << endl;
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
    for (int p = 0; p < 6; ++p) {
        for (int q = 0; q < 6; ++q) {
            for (int r = 0; r < 6; ++r) {
                for (int s= 0; s < 6; ++s) {
                    D = v(p/2,q/2)(r/2,s/2);  // Direct term
                    Ex = v(p/2,q/2)(s/2,r/2); // Exchange term
                    //cout << D << " " << Ex << endl;
                    V(p,q)(r,s) = Z*state(p,q,r,s,D,Ex);
                }
            }
        }
    }

    return V;
}

double HFSolve::state(int p, int q, int r, int s, double D, double Ex){
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

double HFSolve::h0(int alpha,int gamma){
    // the one-boy interaction

    double h = 0;
    if (alpha == gamma){
        double n = alpha/2 + 1.0;
        //cout << "Alpha/n:  " << alpha << "/" << n << endl;
        h = -(Z*Z)/(2*n*n);
        //cout << "this is Z:" << Z << endl;
        //printf(Z)
    }

    return h;
}



mat HFSolve::HF(mat C, field<mat> V){
    /* Sets up the Hartree-Fock matrix
     * using the coefficients given in C
     */

    mat HFmx;
    HFmx.zeros(C.n_rows, C.n_cols);
    for (int alpha = 0; alpha < 6; ++alpha) {
        for (int gamma = 0; gamma < 6; ++gamma) {
            double interaction = 0;
            for (int p = 0; p < N; ++p) {
                for (int beta = 0; beta < 6; ++beta) {
                    for (int delta = 0; delta < 6; ++delta) {
                        interaction += C(p,beta)*C(p,delta)*V(alpha,beta)(gamma,delta);
                        //interaction = interaction + 0*alpha + 0*beta + 0*gamma + 0*delta + 0*p;
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

    //cout << "diff.max() : " << abs(e_v.max() - e_v_prev.max()) << endl;
    while (abs(e_v.min() - e_v_prev.min()) > tolerance){ // convergence test
        iters = iters + 1;
        //cout << "in while iteration: " << iters << endl;
        //for (int i = 0; i < 6; ++i) {
        //    cout << e_v_prev(i) << " " << e_v(i) << endl;
        //}
        e_v_prev = e_v;
        // return the eigenvalues of the HF-mx to e_v and the eigenvectors to C.
        eig_sym(e_v,C,HF(C,V));
        C = trans(C);

    }
    cout << "------------------------------" << endl;
    cout << "iterations: " << iters << endl;
    //cout << "eigenvalues: " << endl;
    int minIndex = 0;
    for (int i = 0; i < 6; ++i) {
        if (e_v[i] < e_v[minIndex]){
            minIndex = i;
        }
        //cout << e_v[i] << endl;
    }
    double E = calc_energy(C,V);
    cout << "Ground State Energy: " << E << endl;


}

double HFSolve::calc_energy(mat C, field<mat> V){

    mat Ci = zeros(6,6);
    double s;
    for(int p= 0; p<6;p++){
        for(int q=0; q<6; q++){
            s = 0;
            for(int k=0; k<N; k++){
                s+=C.at(k,p) * C.at(k,q);
                Ci.at(p,q) = s;
            }
        }
    }
    double Energy = 0;
    for(int alpha = 0; alpha < 6; alpha++){
        Energy += h0(alpha,alpha)*Ci.at(alpha,alpha);
        for(int beta = 0; beta < 6; beta++){
            for(int gamma = 0; gamma<6; gamma++){
                for(int delta = 0; delta <6; delta++){
                    Energy += 0.5 * Ci.at(alpha, gamma) * Ci.at(beta, delta) * V(alpha, beta)(gamma, delta);
                }
            }
        }
    }

    //double Energy = 0;
    //for (int i = 0; i < N; ++i) {
    //    for (int alpha = 0; alpha < 6; ++alpha) {
    //        for (int beta = 0; beta < 6; ++beta) {
    //            double interaction = 0;
    //            for (int j = 0; j < N; ++j) {
    //                for (int gamma = 0; gamma < 6; ++gamma) {
    //                    for (int delta = 0; delta < 6; ++delta) {
    //
    //                        interaction += C(i,alpha)*C(j,beta)*C(i,gamma)*C(j,delta)*V(alpha,beta)(gamma,delta);
    //                    }
    //                }
    //            }
    //
    //            Energy += C(i,alpha)*C(i,beta)*h0(alpha,beta) + 0*interaction;
    //        }
    //
    //    }
    //}
    return Energy;
}




