#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include <string>            // to_string
#include <boysfunction.h>
#include <basis.h>
#include <integrator.h>
#include <hfsolve.h>
#include <returnhermitecoeffs.h>
#include <kineticenergy.h>
#include <setuphermiteintegral.h>
#include <contracted.h>

double pi = 4*atan(1);

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {
    /* Make sure any basis to read from file is located in build folder!
     *
     * Yet to be implemented (as per 13/3/2014):
     * (1) Correct convergence conditions in HFSolve
     * (2) GTOs
     * (3) The Thijssen implementation of spin
     *
    */


/*
    //set number of protons and electrons
    int Z = 4;   // Number of protons
    int N = 4;   // Number of electrons
    int Ns = 4;  // 6 states

    basis Bs(N, 0);              //creating the basis object
    string filename;
    filename = "m_elements_c.dat";
    Bs.read(filename, Z); //reading basis from file
    Bs.set_orthonormal(true);

//    //Solving for N,Z with the provided basis
    HFSolve object (Z,N);
    double E = object.Solve(Bs);

    cout << "Energy of the ground state= " << E << endl;
    */

    double weight = 1.0;
    double a1,b1;
    int i,j,k,l,m,n;
    vec A,B;

    // PrimitiveA:
    a1 = 0.2;
    weight = 1;
    i = 0;
    k = 0;
    m = 1;
    A = {1.2, 2.3, 3.4};

    // PrimitiveB:
    b1 = 0.3;
    weight = 1;
    j = 1;
    l = 0;
    n = 0;
    B = {-1.3,1.4,-2.4};

    Primitive primitiveA(weight,i,k,m,a1,A);
    Primitive primitiveB(weight,j,l,n,b1,B);

    //Core position
    integrator AB (primitiveA, primitiveB);
    vec3 C = {2.3,0.9,3.2};

    AB.setupRtuv(C);
    Primitive contr[2] = {primitiveA, primitiveB};
    contracted BASE (2,contr);

    return 0;
} // End: function output()
