/* Adding a test.cpp-file to test the implementation of the Class HFSolver.cpp
 * HFSolver should be initialized with a inputfile
 */

#include <iostream>
#include <armadillo>
#include <fstream>
#include <hfsolve.h>
#include <string>

using namespace std;
using namespace arma;

int mainf(){

    string filename;
    filename = "m_elements_c.dat";
    // Make sure the filename is located in the build-folder !!!
    // Yet to be implemented (6.march 2014):
    // (1) A more versatile class-structure
    // (2) The Thijssen spin implementation (reducing the size of the matrix)
    // (3) Gaussian Type Orbitals
    // (4) If possible/sensible; a possibility to calculate excited states using excited SDs. (particle-hole states)
    //     In effect this means to permute the C-matrix after each iteration, so that the two one-particle states is not the two lowest lying states.

    int Z = 2; // two protons
    int N = 2; // two electrons
    int Ns = 6; //6 states
    field<mat> V;
    HFSolve object (Z,N,Ns);
    V = object.init(filename);
    object.Solve(V);
    return 0;
}
