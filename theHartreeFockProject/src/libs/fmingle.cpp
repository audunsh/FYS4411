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
    basisbank wrapped (BS);
    fminglebasisbank = wrapped;
}

void fmingle::add_nucleus(vec3 corePos, int charge){
    fminglebasisbank.bs.add_nucleus(corePos, charge);
}

void fmingle::add_orbitals(vec3 corePos, string config){
    try{
        //fminglebasisbank.add_
    }
    catch(int e){
        cout << "Could not find requested basis:" << config << endl;
    }
}

void fmingle::rhf_solve(int nElectrons){

    if(initialized == 0){
        fminglebasisbank.bs.set_size(fminglebasisbank.bs.Nstates);
        fminglebasisbank.bs.init_integrals();
        initialized =1;
    }
    rhfsolve solverobject(fminglebasisbank.bs, nElectrons);
    fminglesolver_rhf =solverobject;
    rhf_energy = fminglesolver_rhf.solve();
    cout << "-------------------------------------------------------------------" << endl;
    cout << "  Restricted Hartree-Fock energy:" << rhf_energy << endl;
    cout << "-------------------------------------------------------------------" << endl;
};

void fmingle::uhf_solve(int nElectronsUp, int nElectronsDown){
    if(initialized == 0){
        fminglebasisbank.bs.set_size(fminglebasisbank.bs.Nstates);
        fminglebasisbank.bs.init_integrals();
        initialized =1;
    }
    uhfsolve solverobject(fminglebasisbank.bs, nElectronsUp, nElectronsDown);
    fminglesolver_uhf =solverobject;
    uhf_energy = fminglesolver_uhf.solve();
    cout << "-------------------------------------------------------------------" << endl;
    cout << "Unrestricted Hartree-Fock energy:" << uhf_energy << endl;
    cout << "-------------------------------------------------------------------" << endl;
};

void fmingle::ccd_solve(int nElectrons){
}

void fmingle::ccsd_solve(int nElectrons){

}
