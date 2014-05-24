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
    if(false){
        //Calculate Beryllium Ground state using STO-3G

        /*
        basis BS;
        int nElectrons = 2;
        double nProtons = 2;
        BS.init_H2({2,1.3,0},{2,2.7,0}); //insert parameter dist here (calculation is however still off for molecules)
        BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction
        hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
        E = object.solve();                          //solve for the given basis
        cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (Approx. " << 27.212*E << " eV)" << endl;
        */

        basis BS;
        int nElectrons = 4;
        double nProtons = 4;
        BS.init_STO_3G("Be", nProtons);
        BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction
        BS.h.print();
        BS.v.print();
        BS.S.print();
        hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
        double E = object.solve();                          //solve for the given basis
        cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (Approx. " << 27.212*E << " eV)" << endl;
    }
    if(true){

        //Calculate H2 Ground state
        basis BS;               //initialize basis object
        int N = 30;  //Grid to calculate is NxN
        mat energies;
        energies.zeros(N,N);

        hartreefocksolver object (BS, 2,1);
        vec3 corePosH1, corePosH2;
        vec3 molecularCenter = {0,0,0};

        double x = 0;
        double x0 = 0.1;
        double dx = 0.1;

        double y = 0;
        double y0 = 0;
        double dy = 0.1;

        vec3 dB1, dB2;

        for(int i=0; i<N;i++){
            for(int j=0; j<N;j++){
                x = i*dx + x0;
                y = j*dy + y0;

                corePosH1 = {0,0,0};
                corePosH2 = {x,y,0};

                dB1 = corePosH1 + molecularCenter;
                dB2 = corePosH2 + molecularCenter;

                BS = basis(); //reinitializing class

                BS.init_H2(dB1,dB2);
                BS.init_integrals();

                object = hartreefocksolver(BS,2,2); //reinitializing class



                energies(i,j) = object.solve();
                cout << "series: [ " << i << " | " << j << " ]  " <<   " Energy converge occurs at at " << energies(i,j) << ". Distance: " << sqrt(x*x+y*y) << endl;
            }
        }
        energies.save("H2_004.dataset", raw_ascii);
        cout << "Calculation complete, file saved to disk." << endl;

    }


    if(false){
        //Perform a lowest eigenenergy fit of a H2Be molecule
        basis BS;               //initialize basis object
        int N = 100;
        double dist;
        mat energies;
        energies.zeros(N,N);
        hartreefocksolver object (BS, 6,6);
        vec3 corePosH1, corePosH2, corePosO;
        vec3 molecularCenter = {0,0,0};

        double halfDist = 0;
        double halfDist0 = 0.5;  //1*sin(0.91);//-0.5;
        //double halfDist0 = 1.00;
        double d_halfDist = 0.02;

        double ODist = 0;
        double ODist0 = 0.0;//1*cos(0.91);//-0.5;
        //double ODist0 = 0.50;
        double d_ODist = 0.02;


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
                corePosO =  {0,0,ODist};

                dB1 = corePosH1;// + molecularCenter;
                dB2 = corePosH2;// + molecularCenter;
                dB3 = corePosO ;// + molecularCenter;
                //cout << "       " << halfDist << " " << ODist << endl;
                //dB1.print();
                //cout << endl;
                //dB2.print();
                //cout << endl;
                //dB3.print();
                //cout << endl;
                cout << "              r:" << sqrt(halfDist*halfDist + ODist*ODist) << endl;


                BS.init_H2Be(dB1,dB2,dB3);
                BS.init_integrals();
                //cout << BS.Nstates << endl;
                //cout << BS.nnInteraction() << endl;
                //cout << BS.evaluateContracted(0+j, {0,0,0}) << endl;
                //cout << BS.evaluateContracted(1+j, {0,0,0}) << endl;
                //cout << BS.evaluateContracted(2+j, {0,0,0}) << endl;


                object.reset(BS,6,6);
                //hartreefocksolver object(BS, 10,8);
                energies(i,j) = object.solve();
                cout << "series:" << i << " -->  At angle:" << setprecision(10) << 2*acos(ODist/sqrt(ODist*ODist + halfDist*halfDist)) << " the energy converges at " << energies(i,j) << endl;
                //cout << "       " << halfDist << " " << ODist << endl;
                cout << "-------------------" << endl;
            }
        }
        energies.print();
        energies.save("H2Be_004.dataset", raw_ascii);
        double E = energies(0,0);
        cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (" << 27.212*E << " eV)" << endl;        //print out approximated ground state energy
    }




    return 0;

} // End: function output()

