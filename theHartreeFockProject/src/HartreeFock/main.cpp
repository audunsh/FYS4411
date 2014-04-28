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
/*
    double weight = 1.0;
    int exponent = 2.0;
    int xExponent = 2.0;
    int yExponent = 2.0;
    int zExponent = 2.0;
    vec AnucleusPosition = {0,0,0};
    vec BnucleusPosition = {1,1,1};
    integrator integral;
    Primitive a(weight,xExponent,yExponent,zExponent,exponent,AnucleusPosition);
    Primitive b(weight,xExponent,yExponent,zExponent,exponent,BnucleusPosition);

    double Sab = integral.overlapIntegral(a,b);
    cout << Sab << endl;
*/
    double weight = 1.0;
    double a1,b1;
    int i,j,k,l,m,n;
    vec A,B;

    // PrimitiveA:
    a1 = 0.2;
    weight = 1;
    i = 2;
    m = 2;
    k = 1;
    A = {1.2, 2.3, 3.4};

    // PrimitiveB:
    b1 = 0.3;
    weight = 1;
    j = 0;
    l = 1;
    n = 2;
    B = {-1.3,1.4,-2.4};

    Primitive primitiveA(weight,i,k,m,a1,A);
    Primitive primitiveB(weight,j,l,n,b1,B);

    // testing ReturnHermiteCoeffs:
    integrator AB (primitiveA, primitiveB);


    /*
    //ReturnHermiteCoeffs Coeffs;

    field <cube> Eab = Coeffs.ReturnCoeffs(primitiveA,primitiveB);

    //double Sab2 = Eab(0)(i,j,0)*Eab(1)(k,l,0)*Eab(2)(m,n,0)*pow(pi/(a1+b1),3.0/2);
    cout << setprecision(13) << " Sab = " << Sab2 << " and should be:" << endl;
    cout << " Sab = -0.07329386373895" << endl;

    // testing kinetic energy
    KineticEnergy T(Eab,&primitiveA,&primitiveB);

    double Tab = T.ReturnKineticIntegral();

    cout << "-----------------------------------------------------------" << endl;
    cout << "Tab= " << Tab << " And should be:" << endl;
    cout << "Tab= -0.01598401092187" << endl;


    // testing the Nuclei-Electron integral
    int t,u,v;
    double p = a1+b1;
    vec corePosition = {0.0, 0.0, 0.0};
    setupHermiteIntegral HermiteIntegral(primitiveA,primitiveB,corePosition);
    field <cube> Rtuv = HermiteIntegral.ReturnHermiteIntegral();
    cout << Rtuv(n)(t+1,u+1,v+1) << endl;

    */

    return 0;
} // End: function output()
