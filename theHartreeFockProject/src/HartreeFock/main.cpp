#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include <string>            // to_string
#include <basis.h>
#include <boysfunction.h>
#include <integrator.h>
#include <hfsolve.h>
#include <returnhermitecoeffs.h>

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


    //set number of protons and electrons
    int Z = 4;   // Number of protons
    int N = 4;   // Number of electrons
    int Ns = 4;  // 6 states

    basis Bs(N, 0);              //creating the basis object
    string filename;
    filename = "/home/goranbs/goran/CompPhys/FYS4411\ -\ CompPhys2/build-theHartreeFockProject-Desktop-Release/testingHF/m_elements_c.dat";
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
    i = m = 0;
    k = 2;
    A = {1.2, 2.3, 3.4};

    // PrimitiveB:
    b1 = 0.3;
    weight = 1;
    j = l = 1;
    n = 0;
    B = {-1.3,1.4,-2.4};

    Primitive primitiveA(weight,i,k,m,a1,A);
    Primitive primitiveB(weight,j,l,n,b1,B);
    // testing integrator:

    //integrator integral2;
    //double Sab2 = integral2.overlapIntegral(primitiveA,primitiveB);
    //cout << " Sab = " << Sab2 << " and should be -7.329386373895e-02" << endl;

    // testing ReturnHermiteCoeffs:


    ReturnHermiteCoeffs Coeffs;

    field <cube> Eab = Coeffs.ReturnCoeffs(primitiveA,primitiveB);

    double Sab2 = Eab(0)(i,j,0)*Eab(1)(k,l,0)*Eab(2)(m,n,0)*pow(pi/(a1+b1),3.0/2);
    cout << setprecision(13) << " Sab = " << Sab2 << " and should be:" << endl;
    cout << " Sab = -0.07329386373895" << endl;



    // Find the Kinetic energy elements:
    //field <mat> T = Coeffs.ReturnKineticMatrix();

    // print out the Kinetic enegy matrix:
/*
    for (int cor = 0; cor < 3; ++cor) {
        cout << "--------------------------------------------" << endl;
        cout << "------------------- T" << cor << " ---------------------" << endl;
        cout << "--------------------------------------------" << endl;
        for (int iA = 0; iA < T(cor).n_rows; ++iA) {
            for (int iB = 0; iB < T(cor).n_cols; ++iB) {
                cout << T(cor)(iA,iB) << " " ;
            }
            cout << endl;
        }
    }
*/

    //double Tab = Coeffs.ReturnKineticIntegral();

    cout << "------------------------------" << endl;
    //cout << "Tab= " << Tab << " And should be:" << endl;
    cout << "Tab= -0.01598401092187" << endl;
    return 0;
} // End: function output()
