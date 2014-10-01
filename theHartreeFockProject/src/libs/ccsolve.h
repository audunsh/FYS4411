#ifndef CCSOLVE_H
#define CCSOLVE_H

#include <hartreefocksolver.h>
#include <armadillo>

class ccsolve
{
public:
    ccsolve();
    ccsolve(hartreefocksolver object, int nElect);

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
    bool unconverged(double tolerance);
    double CCDQ(int a, int b, int i, int j);
    double CCDL(int a, int b, int i, int j);
    double GetCoupledElement(int a, int b, int c, int d);
    double GetUncoupledElement(int a, int b);
    double energy();
    double equalfunc(int a, int b);
    double CCSD_Single(int i, int a);
    double CCSD_Double(int i, int j, int a, int b);

    //Variables
    int nElectrons;
    int nStates;
    int nFermi;

private:
    field<mat> vmin; //The coupled minimized matrix elements
    field<mat> temp_mo; //quarter sized molecular orbital elements
    mat fmin; //The uncoupled minimized matrix elements
    hartreefocksolver hfobject;
    mat Cm;

    //amplitude tensors
    mat t1;
    field<mat> t2;

    mat t1new;
    field<mat> t2new;
    double eprev;

};

#endif // CCSOLVE_H
