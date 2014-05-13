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

    double nProtons  = 2; //number of protons
    int nElectrons= 2;      //number of electrons
    basis BS;               //initialize basis object

    //Enable line below to init hydrogenlike orbit (precomputed)

    //From the first part of the project, we obtained He: -2.8315, Be:-14.5055 using 3 of the orbits from the basis below
    //-2.807, -14.35188 (fra dragly)
    //BS.init_HTO4(nProtons); //set up hydrogenlike basis

    //minor testing below
    /*
    Primitive S1A(0.44063,0,0,0,6.36242139,{0,0,0});
    Primitive S1B(0.426158,0,0,0,1.15892300,{0,0,0});
    vec3 nucleiPos = {0,0,0};
    BoysFunction boys(3);
    integrator AB(S1A,S1A, boys);
    AB.setupRtuv(nucleiPos);
    cout << AB.kinetic() << endl;
    cout << AB.pNuclei() << endl;
    cout << nProtons*AB.pNuclei() << endl;
    */




    //Enable the two lines below for STO-3G:Be basis
    //For Helium using STO_3G; uncoupled: -1.9317, coupled: 1.0557
    BS.init_STO_3G("He", nProtons); //initialize the STO-3G basis for the Beryllium atom
    BS.init_integrals();  //set up and solve the needed integrals to calculate overlapmatrix, single-particle interaction and two-particle interaction

    //BS.turboCorrection();
    hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
    double E = object.solve();                          //solve for the given basis
    cout << "Ground state energy:" << E << endl;        //print out approximated ground state energy

    return 0;
} // End: function output()
