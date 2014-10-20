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
    //Testing UHF procedure
    if(false){
        cout << "Performing an UHF on Li using STO-6G set" << endl;
        vec3 corePos0 = {0,0,0};                            //setting up some position vectors for the cores
        fmingle myparty1;                                   //creating the system
        myparty1.add_nucleus(corePos0, 3);                  //creating a nucleus with charge 1
        myparty1.fminglebasisbank.add_6_311G2df2pd_li(corePos0);  //creating an electrons positioned at core 1 using STO-6G basis set
        myparty1.uhf_solve(1,2);                            //perform a unrestricted hartree fock procedure for 1 up electron, 2 down electrons
    }

    if(true){
        cout << "Performing UHF+CCSD on H2 using STO-6G set" << endl;
        vec3 corePos1 = {0,0,0};                            //setting up some position vectors for the cores
        vec3 corePos2 = {0,0,1.4};
        fmingle myparty2;                                   //creating the system
        myparty2.printing = true;
        myparty2.add_nucleus(corePos1, 1);                  //creating a nucleus with charge 1
        myparty2.add_nucleus(corePos2, 1);                //creating a nucleus with charge 1
        myparty2.fminglebasisbank.add_STO_6G_h(corePos1);  //creating an electrons positioned at core 1 using STO-6G basis set
        myparty2.fminglebasisbank.add_STO_6G_h(corePos2); //creating an electrons positioned at core 2 using STO-6G basis set
        myparty2.uhf_solve(1,1);                            //perform a unrestricted hartree fock procedure for 1 up electron, 2 down electrons
        myparty2.ccsd_solve(2);

    }

    if(false){
        cout << "Performing RHF+CCSD on H2 using STO-6G set" << endl;
        vec3 corePos1 = {0,0,0};                            //setting up some position vectors for the cores
        vec3 corePos2 = {0,0,1.4};
        fmingle myparty2;                                   //creating the system
        myparty2.printing = true;
        myparty2.add_nucleus(corePos1, 1);                  //creating a nucleus with charge 1
        myparty2.add_nucleus(corePos2, 1);                //creating a nucleus with charge 1
        myparty2.fminglebasisbank.add_STO_6G_h(corePos1);  //creating an electrons positioned at core 1 using STO-6G basis set
        myparty2.fminglebasisbank.add_STO_6G_h(corePos2); //creating an electrons positioned at core 2 using STO-6G basis set
        myparty2.rhf_solve(2);                            //perform a unrestricted hartree fock procedure for 1 up electron, 2 down electrons
        myparty2.ccsd_solve(2);

    }

    return 0;

}

