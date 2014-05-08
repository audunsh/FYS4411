#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include <string>
#include <boysfunction.h>
#include <basis.h>
#include <integrator.h>
#include <hfsolve.h>
#include <contracted.h>

double pi = 4*atan(1);

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {
    /* The code seems to be running fine, despite the result being obviously faulty.
     * The error arise in the solver itself (HFSolve) as we still havent implemented
     * the correct use of the density matrix.
     *
     * Pending work is then to work through the solver itself step by step, following Thijessen closely.
     *
     * This should not take too long.
     *
     * Audun, 8/5/14
     */


    /*
    //set number of protons and electrons
    int Z = 4;   // Number of protons
    int N = 4;   // Number of electrons
    int Ns = 4;  // 6 states

    basis BS(N, 0);              //creating the basis object
    string filename;
    filename = "m_elements_c.dat";
    BS.read(filename, Z); //reading basis from file
    BS.set_orthonormal(true);

    cout << "Energy of the ground state= " << E << endl;
    */
    basis BS(3); //set up a basis containing 3 contracted/orbitals
    BS.init_STO_3G("Be"); //initialize the STO-3G basis for the Beryllium atom
    BS.init_integrals();  //set up and solve the needed integrals to calculate overlapmatrix, single-particle interaction and two-particle interaction
    HFSolve object (4,3); //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
    double E = object.Solve(BS); //solve for the given basis
    cout << "Energy of the ground state= " << E << endl; //print out approximated ground state energy
    return 0;
} // End: function output()
