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

        //From the first part of the project, we obtained He: -2.8315, Be:-14.5055 using 3 of the orbits from the basis below
        //-2.807, -14.35188 (fra dragly)

        //Enable the two lines below for STO-3G:Be basis

        //BS.init_STO_3G("He", nProtons); //initialize the STO-3G basis for the Beryllium atom (ion 2+ in current config)

        //BS.init_molecule("O", {8}, {0,0,0});

        //BS.init_Be2({0,0,0},{0,CoreDist,0}); //insert parameter dist here (calculation is however still off for molecules)
        //BS.init_H2({0,0,0},{0,CoreDist,0}); //insert parameter dist here (calculation is however still off for molecules)

        //string filename = "";
        //BS.init_Be2({2,1.3,0},{2,2.7,0}); //insert parameter dist here (calculation is however still off for molecules)
        //BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction
        //hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
        //E = object.solve();                          //solve for the given basis
        //cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (" << 27.212*E << " eV)" << endl;



        int N = 50;
        double dist;
        mat energies;
        energies.zeros(N,N);
        hartreefocksolver object (BS, 10,8);
        vec3 corePosH1, corePosH2, corePosO;
        vec3 molecularCenter = {1,1,0};

        double halfDist = 0;
        double halfDist0 = 0.05;  //1*sin(0.91);//-0.5;
        //double halfDist0 = 1.00;
        double d_halfDist = 0.2;

        double ODist = 0;
        double ODist0 = 0.05;//1*cos(0.91);//-0.5;
        //double ODist0 = 0.50;
        double d_ODist = 0.2;


        vec3 dB1, dB2, dB3;

        for(int i=0; i<N;i++){
            for(int j=0; j<N;j++){

                halfDist = i*d_halfDist + halfDist0;
                ODist =    j*d_ODist    + ODist0;    //halfDist*sin(0.660796);

                //halfDist = j*d_ODist   *sin(0.91)+halfDist0;
                //ODist =    j*d_ODist   *cos(0.91)+ODist0;    //halfDist*sin(0.660796);

                //ODist = (i+1)*d_ODist;
                corePosH1 = {0,-halfDist,0};
                corePosH2 = {0, halfDist,0};
                corePosO =  {ODist,0,0};

                dB1 = corePosH1 + molecularCenter;
                dB2 = corePosH2 + molecularCenter;
                dB3 = corePosO  + molecularCenter;
                cout << "       " << halfDist << " " << ODist << endl;
                //dB1.print();
                //cout << endl;
                //dB2.print();
                //cout << endl;
                //dB3.print();
                //cout << endl;
                cout << "              r:" << sqrt(halfDist*halfDist + ODist*ODist) << endl;


                BS.init_H2O(dB1,dB2,dB3);
                BS.init_integrals();
                //cout << BS.evaluateContracted(0+j, {0,0,0}) << endl;
                //cout << BS.evaluateContracted(1+j, {0,0,0}) << endl;
                //cout << BS.evaluateContracted(2+j, {0,0,0}) << endl;


                object.reset(BS,10,8);
                //hartreefocksolver object(BS, 10,8);
                energies(i,j) = object.solve();
                cout << "series:" << i << " -->  At angle:" << setprecision(10) << 2*acos(ODist/sqrt(ODist*ODist + halfDist*halfDist)) << " the energy converges at " << energies(i,j) << endl;
                //cout << "       " << halfDist << " " << ODist << endl;
                cout << "-------------------" << endl;
            }
        }
        energies.print();
        energies.save("angles_28_001.dataset", raw_ascii);
        double E = energies(0,0);
        cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (" << 27.212*E << " eV)" << endl;        //print out approximated ground state energy



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

