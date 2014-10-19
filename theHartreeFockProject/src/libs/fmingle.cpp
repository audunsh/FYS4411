#include "fmingle.h"
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

fmingle::fmingle()
{
    //this class will be the main user interface for the solver
    //A typical runtime should only communicate with the respective solvers through this class
    basis BS;
    initialized = 0;
    printing = false;
    basisbank wrapped (BS);
    fminglebasisbank = wrapped;
    fminglebasisbank.bs.Nstates = 0;
    report = "";
}

void fmingle::add_nucleus(vec3 corePos, int charge){
    fminglebasisbank.bs.add_nucleus(corePos, charge);
}

void fmingle::add_orbitals(vec3 corePos, string config){
}

void fmingle::initialize(){
    fminglebasisbank.bs.set_size(fminglebasisbank.bs.Nstates);
    fminglebasisbank.bs.init_integrals();
    HFSolve solverobject(fminglebasisbank.bs);
    fminglesolver_hf = solverobject;
}

void fmingle::rhf_solve(int nElectrons){
    if(initialized == 0){initialize();
    }
    initialized =1;
    fminglesolver_hf.solve_rhf(nElectrons);
    rhf_energy = fminglesolver_hf.energy;
    if(printing){
    cout << "-------------------------------------------------------------------" << endl;
    cout << std::setprecision(14)<<"  Restricted Hartree-Fock energy:" << rhf_energy << endl;
    cout << "-------------------------------------------------------------------" << endl;}
};

void fmingle::uhf_solve(int nElectronsUp, int nElectronsDown){
    if(initialized == 0){initialize();
    }
    initialized =2;
    fminglesolver_hf.solve_uhf(nElectronsUp, nElectronsDown);
    uhf_energy = fminglesolver_hf.energy;
    if(printing){
    cout << "-------------------------------------------------------------------" << endl;
    cout << std::setprecision(14)<<"Unrestricted Hartree-Fock energy:" << uhf_energy << endl;
    cout << "-------------------------------------------------------------------" << endl;}
};

void fmingle::ccsd_solve(int nElectrons){
    if(initialized==0){
        cout << "The basis is not initialized." << endl;
        //Do nothing.
    }
    if(initialized==1){
        //Perform ccd for a RHF basis
        ccsolve solverobject(fminglesolver_hf);
        fminglesolver_cc = solverobject;
        fminglesolver_cc.init_RHF_basis();
        fminglesolver_cc.CCSD(nElectrons);
        if(printing){
        cout << "-------------------------------------------------------------------" << endl;
        cout << std::setprecision(14)<< "       CCSD Electron correlation:" << fminglesolver_cc.correlation_energy << endl;
        cout << "-------------------------------------------------------------------" << endl;

        cout << "-------------------------------------------------------------------" << endl;
        cout << std::setprecision(14)<<"                    Total energy:" << fminglesolver_cc.correlation_energy+rhf_energy << endl;
        cout << "-------------------------------------------------------------------" << endl;}
        correlation_energy = fminglesolver_cc.correlation_energy;


    }
    if(initialized==2){
        //Perform ccd for a RHF basis
        //cout << "Entering UHF CCSD procedure."<< endl;
        ccsolve solverobject(fminglesolver_hf);
        fminglesolver_cc = solverobject;
        fminglesolver_cc.init_UHF_basis();
        fminglesolver_cc.CCSD(nElectrons);
        correlation_energy = fminglesolver_cc.correlation_energy;
        if(printing){
        cout << "-------------------------------------------------------------------" << endl;
        cout << std::setprecision(14)<<"       CCSD Electron correlation:" << fminglesolver_cc.correlation_energy << endl;
        cout << "-------------------------------------------------------------------" << endl;
        cout << "-------------------------------------------------------------------" << endl;
        cout << std::setprecision(14)<<"                    Total energy:" << fminglesolver_cc.correlation_energy+uhf_energy << endl;
        cout << "-------------------------------------------------------------------" << endl;}
    }


}

void fmingle::ccd_solve(int nElectrons){

}
