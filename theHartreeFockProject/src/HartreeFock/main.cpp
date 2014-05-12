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
#include <hartreefocksolver.h>

double pi = 4*atan(1);

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {
    /**************************************************************************************************************
    // my thought is to make a python program that calls theHartreeFockProject with values Z,N,Ns and a distance,
    // this distance is the absoulute distance between two cores. We run the python program for different distances,
    // and makes a plot of the potential distribution we get from this :-)

    // such a python program is added in the new folder python_programs.

    int Z,N,Ns;
    double dist;
    //test

    Z = (int) argv[1];
    N = (int) argv[2];
    Ns = (int) argv[3];
    dist = (double) argv[4];
    *****************************************************************************************************************/

    double nProtons  = 4; //number of protons
    int nElectrons= 4;    //number of electrons
    basis BS; //initialize basis object

    BS.init_HTO4(nProtons); //set up hydrogenlike basis

    //Enable the two lines below for STO-3G:Be basis
    //BS.init_STO_3G("Be", nProtons); //initialize the STO-3G basis for the Beryllium atom
    //BS.init_integrals();  //set up and solve the needed integrals to calculate overlapmatrix, single-particle interaction and two-particle interaction

    hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
    double E = object.solve();                          //solve for the given basis
    cout << "Ground state energy:" << E << endl;        //print out approximated ground state energy

    return 0;
} // End: function output()
