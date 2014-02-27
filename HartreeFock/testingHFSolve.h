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
    // Hydrogen atom:
    int Z = 2; // two protons
    int N = 2; // two electrons
    int Ns = 6; //6 states
    field<mat> V;
    HFSolve object (Z,N,Ns);
    V = object.init(filename);
    object.Solve(V);
    return 0;
}
