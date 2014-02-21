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
    // Hydrogen atom:
    int Z = 2; // two protons
    int N = 2; // two electrons
    field<mat> V;
    HFSolve object;
    V = object.init(filename,Z,N);

    cout << "---------------------------" << endl;
    cout << "Here comes V" << endl;
    for (int p = 0; p < 6; ++p) {
        for (int q = 0; q < 6; ++q) {
            for (int r = 0; r < 6; ++r) {
                for (int s = 0; s < 6; ++s) {
                    cout << p << " " << q << " " << r  << " " << s << " " << V(p,q)(r,s) << endl;
                }
            }
        }
    }
    cout << "---------------------------" << endl;

    // now we have created V. Then it is time to call the solver

    object.Solve(V);


    return 0;
}
