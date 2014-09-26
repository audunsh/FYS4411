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
    void initT2();
    double GetCoupledElement(int a, int b, int c, int d);
    double GetUncoupledElement(int a, int b);
    //Variables
    int nElectrons;
    int nStates;
    int nFermi;

private:
    field<mat> vmin; //The coupled minimized matrix elements
    mat fmin; //The uncoupled minimized matrix elements
    hartreefocksolver hfobject;

    //amplitude tensors
    mat t1;
    field<mat> t2;

};

#endif // CCSOLVE_H
