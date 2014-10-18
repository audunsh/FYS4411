#ifndef FMINGLE_H
#define FMINGLE_H
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

class fmingle
{
public:
    fmingle();
    void add_nucleus(vec3 corePos, int charge);
    void add_orbitals(vec3 fPos, string config);
    void rhf_solve(int nElectrons);
    void uhf_solve(int nElectronsUp, int nElectronsDown);
    void ccd_solve(int nElectrons);
    void ccsd_solve(int nElectrons);
    int initialized;
    //basis fminglebasis;
    basisbank fminglebasisbank;
    HFSolve fminglesolver;
    rhfsolve fminglesolver_rhf;
    uhfsolve fminglesolver_uhf;
    ccsolve fminglesolver_cc;
    double rhf_energy, uhf_energy, correlation_energy;


};

#endif // FMINGLE_H