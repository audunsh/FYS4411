#ifndef CCSOLVE_H
#define CCSOLVE_H

#include <rhfsolve.h>
#include <hfsolve.h>
#include <armadillo>

class ccsolve
{
public:
    ccsolve();
    ccsolve(HFSolve object, int nElect);

    //Functions
    void SetupMinimizedBasis();
    void SetupT1();
    void SetupT2();
    void CCD();
    void CCSD();
    void initT1();
    void initT2();
    void expandC();
    void ExpandMinimizedBasis();
    void ReportEnergy();
    bool unconverged(double tolerance);
    double CCDQ(int a, int b, int i, int j);
    double CCDL(int a, int b, int i, int j);

    double CCDQ2(int a, int b, int i, int j, field<mat> tf, double Qac, double Qbc, double Qcc, double Qdc);
    double CCDL2(int a, int b, int i, int j, field<mat> tf, double L1ac, double L1bc,double L2ac,double L2bc,double L2cc);
    void retranslate();

    double GetCoupledElement(int a, int b, int c, int d);
    double GetUncoupledElement(int a, int b);
    double energy(field<mat> tf, mat t1f);
    double CCDenergy(field<mat> tf);
    double equalfunc(int a, int b);
    double CCSD_Single(int a, int i);
    double CCSD_Double(int a, int b,int i, int j);

    //Variables
    int nElectrons;
    int nStates;
    int nFermi;
    int counter;
    int maxiter;

    double correlation_energy;

private:
    field<mat> vmin; //The coupled minimized matrix elements
    field<mat> temp_mo; //quarter sized molecular orbital elements
    mat fmin; //The uncoupled minimized matrix elements
    HFSolve hfobject;
    mat Cm;

    //amplitude tensors
    mat t1;
    field<mat> t2;

    mat t1new;
    mat t1c; //t1 current
    mat t1p; //t1 previous


    field<mat> tf;
    field<mat> t1f;

    field<mat> t2new;

    field<mat> t2newL;
    field<mat> t2newQ;
    field<mat> t2new0;

    field<mat> t2newL2a;
    field<mat> t2newL2b;
    field<mat> t2newL2c;

    field<mat> t2newQa;
    field<mat> t2newQb;
    field<mat> t2newQc;
    field<mat> t2newQd;

    field<mat> t2newLp;
    field<mat> t2newQp;
    field<mat> t2new0p;

    field<mat> t2newL2ap;
    field<mat> t2newL2bp;
    field<mat> t2newL2cp;

    field<mat> t2newQap;
    field<mat> t2newQbp;
    field<mat> t2newQcp;
    field<mat> t2newQdp;

    field<mat> t20;
    field<mat> t2c; //t2 current
    field<mat> t2p; //t2 previous


    double eprev;

};

#endif // CCSOLVE_H


