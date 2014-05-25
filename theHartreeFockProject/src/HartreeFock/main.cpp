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
     * Hartree-Fock solver with support for gaussian basis sets.
     * Written in C++ by Audun Skau Hansen & Goran Brekke Svaland
     * Spring, 2014
     *
     * Library Armadillo is required for compilation.
     **************************************************************************************************************/
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

    //Some sample calculations
    if(false){
        basis BS;
        int nElectrons = 4;
        double nProtons = 4;
        BS.init_STO_3G("Be", nProtons);
        BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction
        hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
        double E = object.solve();                          //solve for the given basis
        cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (Approx. " << 27.212*E << " eV)" << endl;
    }
    if(true){
        //Calculate O2 Ground state
        basis BS;               //initialize basis object
        int N = 33;  //Grid to calculate is NxN
        mat energies;
        energies.zeros(N,N);
        int nElectrons = 8;
        double nProtons = 8;
        hartreefocksolver object (BS, nElectrons,nProtons);
        vec3 corePosH1, corePosH2;
        vec3 molecularCenter = {0,0,0};

        double x = 0;
        double x0 = 0;
        double dx = 0.1;

        double y = 0;
        double y0 = 0;
        double dy = 0.1;

        vec3 dB1, dB2;

        for(int i=0; i<N;i++){
            for(int j=0; j<N;j++){
                x = i*dx + x0;
                y = j*dy + y0;

                corePosH1 = {1.525,1.525,0};
                corePosH2 = {x,y,0};

                dB1 = corePosH1 + molecularCenter;
                dB2 = corePosH2 + molecularCenter;

                BS = basis(); //reinitializing class

                BS.init_O2(dB1,dB2);
                BS.init_integrals();

                object = hartreefocksolver(BS,nElectrons,nProtons); //reinitializing class



                energies(i,j) = object.solve();
                cout << "series: [ " << i << " | " << j << " ]  " <<   " Energy convergence occurs at " << energies(i,j) << " (a.u.). Distance: " << sqrt(x*x+y*y) << endl;
            }
        }
        energies.save("H2_006.dataset", raw_ascii);
        cout << "Calculation complete, file saved to disk." << endl;

    }


    if(false){
        //Perform a lowest eigenenergy fit of a H2Be molecule
        basis BS;               //initialize basis object
        int N = 40;
        mat energies;
        energies.zeros(N,N);
        hartreefocksolver object (BS, 6,6);
        vec3 corePosH1, corePosH2, corePosO;
        vec3 molecularCenter = {0,0,0};


        double x = 0;
        double x0 = 0.1;
        double dx = 0.1;

        double y = 0;
        double y0 = 0.1;
        double dy = 0.1;

        vec3 dB1, dB2, dB3;

        for(int i=0; i<N;i++){
            for(int j=0; j<N;j++){
                x = i*dx + x0;
                y = j*dy + y0;

                corePosH1 = {0.05,0.05,0};
                corePosH2 = {0.05,3.05,0};
                corePosO =  {x,y,0};

                dB1 = corePosH1 + molecularCenter;
                dB2 = corePosH2 + molecularCenter;
                dB3 = corePosO  + molecularCenter;

                BS = basis();
                BS.init_H2O(dB1,dB2,dB3);
                BS.init_integrals();

                //object.reset(BS,6,6);
                object = hartreefocksolver(BS,10,10); //reinitializing class
                energies(i,j) = object.solve();
                cout << "series: [ " << i << " | " << j << " ]  " << " At angle:" << setprecision(10) << 2*acos(y/sqrt(y*y + x*x)) << " the energy converges at " << energies(i,j) << " at a absolute distance r = " << sqrt(x*x+y*y) << ". (x,y) = " << "(" << x << "," << y << ")" << endl;
            }
            cout << " " << endl;
        }
        //energies.print();
        energies.save("H2_O_007_coreCharge0.dataset", raw_ascii);
        cout << "Calculation complete, file saved to disk." << endl;
        //double E = energies(0,0);
        //cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (" << 27.212*E << " eV)" << endl;        //print out approximated ground state energy
    }




    return 0;

} // End: function output()

