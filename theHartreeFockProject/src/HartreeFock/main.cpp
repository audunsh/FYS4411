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
    //NOTE: Had some trouble with git (messed up merge, sorry) so I some changes might have been revertet.
    //Just update it as you like, only make sure the changes in createDensityMap remains the same as in this commit.

    // argc = lenght of argv[];
    // argv = [__name__, arg1,arg2,arg3,...arg_argc]
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
    // predefined values, can be set using terminal: ./HartreeFock nProtons nElectrons CoreDist
    double nProtons  = 2;    // number of protons
    int nElectrons   = 2;    // number of electrons
    double CoreDist  = 1.4;  // Distance between particles 1 and 2.
    double E = 0.0;
    if (argc == 4) {
        nProtons = atoi(argv[1]);
        nElectrons = atoi(argv[2]);
        CoreDist = atof(argv[3]);
        //cout << "nProtons = " << nProtons << " nElectrons= " << nElectrons << " CoreDist= " << CoreDist << endl;
        basis BS;               //initialize basis object

        BS.init_H2({0,0,0},{CoreDist,0,0}); //insert parameter dist here (calculation is however still off for molecules)
        BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction
        hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals

        E = object.solve();                          //solve for the given basis
        cout << E << endl;
    }
    else{
        basis BS;               //initialize basis object
        //BS.add_atom_STO3G();

        //Enable line below to init hydrogenlike orbit (precomputed), remember to disable gaussian orbits in line 49-
        //BS.init_HTO4(nProtons); //set up hydrogenlike basis

        //todo, basis:
        //(1) Contracted needs to be able to evaluated at given coordinates
        //(2) Set up monte carlo integration in a cubic space

        //From the first part of the project, we obtained He: -2.8315, Be:-14.5055 using 3 of the orbits from the basis below
        //-2.807, -14.35188 (fra dragly)

        //Enable the two lines below for STO-3G:Be basis

        //BS.init_STO_3G("He", nProtons); //initialize the STO-3G basis for the Beryllium atom (ion 2+ in current config)

        //BS.init_molecule("O", {8}, {0,0,0});

        //BS.init_Be2({0,0,0},{0,CoreDist,0}); //insert parameter dist here (calculation is however still off for molecules)
        //BS.init_H2({0,0,0},{0,CoreDist,0}); //insert parameter dist here (calculation is however still off for molecules)

        //string filename = "";
        BS.init_Be2({2,1.3,0},{2,2.7,0}); //insert parameter dist here (calculation is however still off for molecules)
        BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction
        hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
        E = object.solve();                          //solve for the given basis
        cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (" << 27.212*E << " eV)" << endl;
        /*
        int N = 1;
        double dist;
        mat energies;
        energies.zeros(N,N);
        hartreefocksolver object (BS, 10,8);
        vec3 corePosH1, corePosH2, corePosO;

        for(int i=0; i<N;i++){
            for(int j=0; j<N;j++){
                corePosH1 = {0,0,0};
                corePosH1 = {0,0,0};
                corePosO = {0,0,0};
                BS.init_H2O(corePosH1, corePosH2, corePosO);
                BS.init_integrals();
                hartreefocksolver object(BS, 10,8);
                energies(i,j) = object.solve();
            }
        }
        */

        /*
        for(int i=0; i<10;i++){

            BS.init_Be2({2,1.3,0},{2,2.7,0}); //insert parameter dist here (calculation is however still off for molecules)
            BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction

            hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals

            E = object.solve();                          //solve for the given basis
            cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (" << 27.212*E << " eV)" << endl;        //print out approximated ground state energy
            string filename = "Be2_STO3G.dataset";
            //object.createDensityMap(filename);
            filename = "Be2_STO3G_";
            filename += to_string(1);
            filename += ".dataset";
        */

    }

    return 0;

} // End: function output()

