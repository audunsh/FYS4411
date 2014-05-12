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
    /* Update
     *
     * As the old solver had become cryptic after weeks unatteded,
     * I decided to write a new one.
     *
     * The new solver is written in a readable way following Thijssen very closely.
     *
     * It however does not return the correct energy, and it is unclear where the
     * error arises. Using hydrogen-like orbitals and solving for Be yields an energy
     * below the groundstate (-18.2395), obviously exposing some kind of faulty algorithm.
     *
     * Pedning work is therefore:
     * - Fix hartreefocksolve.cpp
     * - Set up tests for all integrals (they have already been informally tested and all of them passed)
     * - When hartreefocksolver is fixed; run it for STO-3G and compare results to hydrogen-likes.
     *
     * Audun,  11/5/14
     *
     */

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

    /*
    //minor tests going on here (particle-particle);
    vec3 corePos1 = { 0.55,1,3};
    vec3 corePos2 = {-0.52,5,6};
    vec3 corePos3 = {-0.53,1,2};
    vec3 corePos4 = { 0.45,2,4};
    double a = 13.0077;
    double b = 0.121949;
    double c = 0.444529;
    double d = 10.0077;

    Primitive A( 1,0,0,0, a,corePos1);
    Primitive B( 1,0,0,0, b,corePos1);
    Primitive C( 1,0,1,0, c,corePos1);
    Primitive D( 1,0,0,0, d,corePos1);


    BoysFunction boys(3);
    integrator AB(A, B, boys);
    cout << AB.pp(C,D)<< endl;
    */

    double nProtons  = 4; //number of protons
    int nElectrons= 4; //number of electrons

    basis BS; //initialize basis object
    //BS.init_HTO4(nProtons); //setting up hydrogenlike basis
    //BS.h.print();


    BS.init_STO_3G("Be", nProtons); //initialize the STO-3G basis for the Beryllium atom
    //BS.init_Be2({1,0,0}, {0,0,0});
    BS.init_integrals();  //set up and solve the needed integrals to calculate overlapmatrix, single-particle interaction and two-particle interaction

    hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
    double E = object.solve();                          //solve for the given basis
    cout << "Ground state energy:" << E << endl; //print out approximated ground state energy

    return 0;
} // End: function output()
