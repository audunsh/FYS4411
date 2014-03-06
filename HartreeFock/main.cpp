#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include <string>      // to_string
#include <lib.h>       //tqli, free_matrix
#include <testingHFSolve.h>


using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {

    string filename;
    filename = "m_elements_c.dat";
    // Make sure the filename is located in the build-folder !!!
    // Yet to be implemented (6.march 2014):
    // (1) A more versatile class-structure
    // (2) The Thijssen spin implementation (reducing the size of the matrix)
    // (3) Gaussian Type Orbitals
    // (4) If possible/sensible; a possibility to calculate excited states using excited SDs. (particle-hole states)
    //     In effect this means to permute the C-matrix after each iteration, so that the two one-particle states is not the two lowest lying states.

    int Z = 4; // two protons
    int N = 4; // two electrons
    int Ns = 6; //6 states
    field<mat> V;
    HFSolve object (Z,N,Ns);
    V = object.init(filename);
    object.Solve(V);

    return 0;
} // End: function output()







