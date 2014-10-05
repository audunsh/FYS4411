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
#include <ccsolve.h>

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
    string atomType;         // Example of atomType: Helium.
    int val = 1;
    double E = 0.0;
    //cout << argc << endl;
    //cout << argv[0] << argv[1] << argv[2] << argv[3] << argv[4] << endl;
    if (argc == 5) {
        nProtons = atoi(argv[1]);
        nElectrons = atoi(argv[2]);
        CoreDist = atof(argv[3]);
        atomType = (argv[4]);
        //cout << "nProtons = " << nProtons << " nElectrons= " << nElectrons << " CoreDist= " << CoreDist << endl;
        basis BS;               //initialize basis object
        //cout << "atomType = " << atomType << endl;

        if (atomType == "H"){
           BS.init_H2({0,0,0},{CoreDist,0,0}); //insert parameter dist here (calculation is however still off for molecules)
        }
        else if (atomType == "He") {
            BS.init_He2({0,0,0},{CoreDist,0,0});
        }
        else if (atomType == "Be"){
            BS.init_Be2({0,0,0},{CoreDist,0,0});
        }
        else if (atomType == "O"){
            BS.init_O2({0,0,0},{CoreDist,0,0});
        }
        else{
            cout << "We have currently not " << atomType << " in our system" << endl;
            cout << "sys.exit..." << endl;
            val = 0;
        }

        if (val == 1){
            BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction
            hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals

            E = object.solve();                          //solve for the given basis
            cout << E << endl;
            //cout << atomType << endl;
        }
    }
    else{
        cout << "No system arguments, performing predefined computations." << endl;
    }
    //Some sample calculations
    if(false){
        //some simple integral tests
        Primitive primitiveA(1.67545, 1, 0, 0, 5.03315,{0,0,0});
        Primitive primitiveB(1.05357, 1, 0, 0, 1.1696,{0,0,0});
        Primitive primitiveC(0.166903, 1, 0, 0, 0.380389,{0,0,0});
        //-1.9610 (Auduns feil)
        //-6.4507 (Milads rett)

        BoysFunction boys(1);
        integrator AB1(primitiveA, primitiveB, boys);
        integrator AB2(primitiveB, primitiveA, boys);
        integrator AB3(primitiveA, primitiveC, boys);
        integrator AB4(primitiveC, primitiveA, boys);
        integrator AB5(primitiveB, primitiveC, boys);
        integrator AB6(primitiveC, primitiveB, boys);
        integrator AB7(primitiveB, primitiveB, boys);
        integrator AB8(primitiveA, primitiveA, boys);
        integrator AB9(primitiveC, primitiveC, boys);

        vec3 corePos = {0,0,0};
        AB1.setupRtuv(corePos);
        AB2.setupRtuv(corePos);
        AB3.setupRtuv(corePos);
        AB4.setupRtuv(corePos);
        AB5.setupRtuv(corePos);
        AB6.setupRtuv(corePos);
        AB7.setupRtuv(corePos);
        AB8.setupRtuv(corePos);
        AB9.setupRtuv(corePos);
        cout << -8*AB1.pNuclei()+AB1.kinetic() << endl;
        cout << -8*AB2.pNuclei()+AB2.kinetic() << endl;
        cout << -8*AB3.pNuclei()+AB3.kinetic() << endl;
        cout << -8*AB4.pNuclei()+AB4.kinetic() << endl;
        cout << -8*AB5.pNuclei()+AB5.kinetic() << endl;
        cout << -8*AB6.pNuclei()+AB6.kinetic() << endl;
        cout << -8*AB7.pNuclei()+AB7.kinetic() << endl;
        cout << -8*AB8.pNuclei()+AB8.kinetic() << endl;
        cout << -8*AB9.pNuclei()+AB9.kinetic() << endl;
        //double Nucl = -8*(AB1.pNuclei()+AB2.pNuclei()+AB3.pNuclei()+AB4.pNuclei()+AB5.pNuclei()+AB6.pNuclei()+AB7.pNuclei()+AB8.pNuclei()+AB9.pNuclei());
        //double Kine = AB1.kinetic() + AB2.kinetic() + AB3.kinetic()+AB4.kinetic()+AB5.kinetic()+AB6.kinetic()+AB7.kinetic()+AB8.kinetic()+AB9.kinetic();
        double Nucl = -8*(AB1.pNuclei()+AB2.pNuclei()+AB3.pNuclei()+AB4.pNuclei()+AB5.pNuclei()+AB6.pNuclei()+AB7.pNuclei()+AB8.pNuclei()+AB9.pNuclei());
        double Kine = AB1.kinetic() + AB2.kinetic() + AB3.kinetic()+AB4.kinetic()+AB5.kinetic()+AB6.kinetic()+AB7.kinetic()+AB8.kinetic()+AB9.kinetic();
        cout << Nucl << " " << Kine << " " << Nucl+Kine << endl;

    }

    if(false){
        basis BS;
        int nElectrons = 8;
        double nProtons =8;
        BS.init_O({0,0,0});
        BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction
        //BS.init_HTO4(nProtons);
        //BS.printAllContracted();
        hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals

        double E = object.solve();                          //solve for the given basis
        cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (Approx. " << 27.212*E << " eV)" << endl;
    }

    if(false){
        //create density plot of H2O
        double xCenter = 6.0;
        double yCenter = 3.2;

        basis BS;
        int nElectrons = 10;
        double nProtons =10;
        BS.init_H2O({xCenter-1.7,yCenter,0},{xCenter+1.7,yCenter,0},{xCenter,yCenter+0.95,0});
        BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction
        //BS.init_HTO4(nProtons);
        //BS.printAllContracted();
        hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals

        double E = object.solve();                          //solve for the given basis
        object.createDensityMap("H2O_density_map4.dataset");
        cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (Approx. " << 27.212*E << " eV)" << endl;
    }

    if(false){
        //PErform a sweep of distances in the H2 config using CCSD
        basis BS;
        int nElectrons = 2;
        int nProtons = 2;

        int nMoves = 80; //number of moves in the sweep
        vec RHFe, CCSDe;
        RHFe.set_size(nMoves);
        CCSDe.set_size(nMoves);
        double dx = 0.2; //steplength
        for(int i = 0; i<nMoves; i++){
            BS = basis();
            vec3 corePosH1 = {0,0,0};
            vec3 corePosH2 = {0,0,(i+1)*dx};
            BS.init_H2(corePosH1, corePosH2);
            //BS.nucleusCharges(0) = 1;
            //BS.nucleusCharges(1) = 1;
            BS.init_integrals();
            hartreefocksolver object (BS, nElectrons, nProtons);
            RHFe(i) = object.solve();
            //cout << "Hartree-Fock energy:" << energy << endl;
            ccsolve ccobject (object, nElectrons);
            CCSDe(i) = ccobject.energy();
        }
        RHFe.print();
        cout << endl;

        CCSDe.print();


        cout << "End of program." << endl;
    }

    if(true){

        //Calculate H2O Ground state energy
        basis BS;
        int nElectrons =10;
        int nProtons = 10;
        double xCenter = 0.1; //0.1 for same coefficient matrix
        double yCenter =0.3;  //0.3 for same coefficient matrix
        BS.init_H2O({xCenter-1.4756110550780858,yCenter+1.079252144093028,0},{xCenter+1.4756110550780858,yCenter+1.079252144093028,0},{xCenter,yCenter,0});

        //BS.init_O2(corePosH1, corePosH2);
        //BS.nucleusCharges(0) = 1;
        //BS.nucleusCharges(1) = 1;
        //BS.nucleusCharges.print();
        BS.init_integrals();
        //BS.h.print();
        hartreefocksolver object (BS, nElectrons, nProtons);
        double energy = object.solve();
        cout << "Hartree-Fock energy:" << energy << endl;
        ccsolve ccobject (object, nElectrons);
        cout << "End of program." << endl;
    }

    if(false){
        //Calculate H2 Ground state energy
        basis BS;
        int nElectrons =2;
        int nProtons = 2;


        vec3 corePosH1 = {0,0,0};
        vec3 corePosH2 = {0,0,1.4};
        BS.init_H2(corePosH1, corePosH2);

        //BS.init_O2(corePosH1, corePosH2);
        //BS.nucleusCharges(0) = 1;
        //BS.nucleusCharges(1) = 1;
        //BS.nucleusCharges.print();
        BS.init_integrals();
        //BS.h.print();
        hartreefocksolver object (BS, nElectrons, nProtons);
        double energy = object.solve();
        cout << "Hartree-Fock energy:" << energy << endl;
        ccsolve ccobject (object, nElectrons);
        cout << "End of program." << endl;
    }
    if(false){
        //Calculate O2 Ground state
        basis BS;               //initialize basis object
        int N = 33;  //Grid to calculate is NxN
        mat energies;
        energies.zeros(N,N);
        int nElectrons = 16;
        double nProtons = 16;
        hartreefocksolver object (BS, nElectrons,nProtons);
        vec3 corePosH1, corePosH2;
        vec3 molecularCenter = {0,0,0};

        double x = 0;
        double x0 = 0;
        double dx = 0.1;

        double y = 0;
        double y0 = 1;
        double dy = 0.3;

        vec3 dB1, dB2;

        for(int i=0; i<N;i++){
            for(int j=0; j<N;j++){
                x = i*dx + x0;
                y = j*dy + y0;

                corePosH1 = {0,.01,0};
                corePosH2 = {x,y,0};

                dB1 = corePosH1 + molecularCenter;
                dB2 = corePosH2 + molecularCenter;

                BS = basis(); //reinitializing class

                BS.init_O2(dB1, dB2);
                BS.init_integrals();


                object = hartreefocksolver(BS,nElectrons,nProtons); //reinitializing class



                energies(i,j) = object.solve();
                cout << "series: [ " << i << " | " << j << " ]  " <<   " Energy convergence occurs at " << energies(i,j) << " (a.u.). Distance: " << sqrt(x*x+y*y) << endl;
            }
        }
        energies.save("Be2_100.dataset", raw_ascii);
        cout << "Calculation complete, file saved to disk." << endl;

    }

    if(false){
        //Perform a sweep along a given axis, save array to disk

        basis BS;               //initialize basis object
        int N = 2;             //Number of calculations
        vec energies;
        energies.zeros(N);
        int nElectrons = 4;
        double nProtons = 4;

        hartreefocksolver object (BS, nElectrons,nProtons);
        vec3 corePosH1, corePosH2;
        vec3 molecularCenter = {0,0,0};

        double x = 0;
        double x0 = 5.5;
        double dx = 0.1;
        vec3 dB1, dB2;

        for(int i=0; i<N;i++){
            x = i*dx + x0;

            corePosH1 = {0,0,0};
            corePosH2 = {x,0,0};

            dB1 = corePosH1 + molecularCenter;
            dB2 = corePosH2 + molecularCenter;

            BS = basis(); //reinitializing class

            BS.init_He2(dB1, dB2);
            BS.init_integrals();

            object = hartreefocksolver(BS,nElectrons,nProtons); //reinitializing class

            energies(i) = object.solve();
            cout << "series: [ " << i << " ]  " <<   " Energy convergence occurs at " << energies(i) << " (a.u.). Distance: " << x << endl;

        }
        //energies.save("Sweep_He2_101.dataset", raw_ascii);
        cout << "Calculation complete, file saved to disk." << endl;


    }


    if(false){
        //Perform a lowest eigenenergy fit of a H2O molecule
        basis BS;               //initialize basis object
        int N = 25;
        mat energies;

        int nElectrons = 10;
        double nProtons = 10;
        energies.zeros(N,N);
        hartreefocksolver object (BS, nElectrons,nProtons);
        vec3 corePosH1, corePosH2, corePosO;
        vec3 molecularCenter = {0,0,0};

        double x = 0;
        double x0 = 1.3;
        double dx = 0.02;

        double y = 0;
        double y0 = 0.8;
        double dy = 0.02;

        vec3 dB1, dB2, dB3;

        for(int i=0; i<N;i++){
            for(int j=0; j<N;j++){
                x = i*dx + x0;
                y = j*dy + y0;

                corePosH1 = {-x,0,0};
                corePosH2 = {x,0,0};
                corePosO =  {0,y,0};

                dB1 = corePosH1 + molecularCenter;
                dB2 = corePosH2 + molecularCenter;
                dB3 = corePosO  + molecularCenter;

                BS = basis();
                BS.init_H2O(dB1,dB2,dB3);
                BS.init_integrals();

                //object.reset(BS,6,6);
                object = hartreefocksolver(BS,nElectrons,nProtons); //reinitializing class
                energies(i,j) = object.solve();
                cout << "series: [ " << i << " | " << j << " ]  " << " At angle:" << setprecision(10) << 2*acos(y/sqrt(y*y + x*x)) << " the energy converges at " << energies(i,j) << " at a absolute distance r = " << sqrt(x*x+y*y) << ". (x,y) = " << "(" << x << "," << y << ")" << endl;
            }
            cout << " " << endl;
        }
        //energies.print();
        energies.save("H2O_208.dataset", raw_ascii);
        cout << "Calculation complete, file saved to disk." << endl;
        //double E = energies(0,0);
        //cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (" << 27.212*E << " eV)" << endl;        //print out approximated ground state energy
    }




    return 0;

} // End: function output()

