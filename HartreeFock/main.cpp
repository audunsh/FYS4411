#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include <string>            // to_string
#include <lib.h>             //tqli, free_matrix
#include <basis.h>
#include <boysfunction.h>
#include <integrator.h>
#include <hfsolve.h>


using namespace std;
//using namespace arma;

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
    int Z = 2;   // Number of protons
    int N = 2;   // Number of electrons
    int Ns = 4;  // 6 states

    basis Bs(3, 0);              //creating the basis object
    string filename;
    filename = "m_elements_c.dat";
    Bs.read("m_elements_c.dat", Z); //reading basis from file
    Bs.set_orthonormal(true);

    //Solving for N,Z with the provided basis
    HFSolve object (Z,N);
    object.Solve(Bs);



    return 0;
} // End: function output()







