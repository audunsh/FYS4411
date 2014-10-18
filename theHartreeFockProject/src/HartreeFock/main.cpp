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
#include <rhfsolve.h>
#include <uhfsolve.h>
#include <ccsolve.h>
#include <basisbank.h>
#include <fmingle.h>

double pi = 4*atan(1);


using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {
    /**************************************************************************************************************
     * Fermion Mingle Quantum Solver
     * A solver for quantum many-body problems with support for:
     * - Restricted and unrestricted Hartree Fock
     * - Coupled Cluster (CCD, CCSD, CCSDT*)
     * - Gaussian Basis sets
     * - Fermion density evaluation
     *
     * (*) remains to be implemented
     * Written as a university project in C++ by Audun Skau Hansen & Goran Brekke Svaland | Comp-Phys | UiO | 2014
     *
     * Library Armadillo is required for compilation.
     **************************************************************************************************************/
    vec3 corePos1 = {0,0,0};                            //setting up some position vectors for the cores
    vec3 corePos2 = {0,0,1.4};
    fmingle mysystem;                                   //creating the system
    mysystem.add_nucleus(corePos1, 3);                  //creating a nucleus with charge 1
    //mysystem.add_nucleus(corePos2, 1);                //creating a nucleus with charge 1
    mysystem.fminglebasisbank.add_STO_6G_li(corePos1);  //creating an electrons positioned at core 1 using STO-6G basis set
    //mysystem.fminglebasisbank.add_STO_6G_h(corePos2); //creating an electrons positioned at core 2 using STO-6G basis set
    mysystem.uhf_solve(2,1);                            //perform a unrestricted hartree fock procedure for 1 up electron, 2 down electrons
    mysystem.ccsd_solve(2);                             //perform a coupled cluster procedure for 2 electrons
    cout << mysystem.report;
    return 0;

}

