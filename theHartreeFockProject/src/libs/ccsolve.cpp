#include "ccsolve.h"
#include <hartreefocksolver.h>

ccsolve::ccsolve()
{
}

ccsolve::ccsolve(hartreefocksolver object, int nElect)
{
    /*
     * This class should be initialized with a hartreefocksolver object containing an energy-minimized basis for the Coupled-Cluster (CC) solver.
     * The CC Solver will then perform the following operations
     * 1. Set up a new basis from atomic to molecular orbitals (Done, but do I need also the f_i^a elements
     * 2. Initialize the cluster amplitudes with a given initial value (Begin with CCD only)
     * 3. Solve the amplitude equations
     * 4. Return the CC-energy
    */
    hfobject = object;
    nElectrons = nElect; //fermi level is defined by number of electrons
    nStates = hfobject.C.n_cols;
    cout << "CCSolve initialized." << endl;
    cout << "Identifying number of orbitals in the basis" << endl;
    cout << "Found " << nStates << " number of basis functions in the HF Rothaan expansion." << endl;
    //SetupMinimizedBasis();
    SetupT1();
    SetupT2();
}

void ccsolve::initT2(){
    for(int a=nElectrons; a<nStates; a++){
        for(int b=nElectrons; b<nStates; b++){
            for(int i=0;i<nElectrons;i++){
                for(int j=0;j<nElectrons;j++){
                    t2(a,b)(i,j) = vmin(a,b)(i,j)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                }
            }
        }
    }
}

double ccsolve::CCDQ(){
    double SM = 0;
    double Qa = 0;
    double Qb = 0;
    int a,b,i,j; //These should be included in the function call
    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; d++){
                for(int d = nElectrons; d<nStates; d++){
                    Qa += vmin(k,l)(c,d)*t(c,d)(i,j)*t(a,b)(k,l);
                }
            }
        }
    }
    Qa/=4.0;

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; d++){
                for(int d = nElectrons; d<nStates; d++){
                    Qb += vmin(k,l)(c,d)*(t(a,c)(i,k)*t(b,d)(j,l) - t(a,c)(j,k)*t(b,d)(i,l));
                }
            }
        }
    }



}

void ccsolve::CCD(){
    //Set up and solve the CCD-equation
    initT2(); //Set up initial guess following S-B, p.289
}

void ccsolve::SetupT1(){
    t1.set_size(nStates, nStates);
    t1.zeros();
    cout << "Successfully initialized t1 amplitudes." << endl;
}

void ccsolve::SetupT2(){
    t2.set_size(nStates, nStates);
    for(int i = 0; i<nStates; i++){
        for(int j =0;j<nStates; j++){
            t2(i,j).set_size(nStates, nStates);
            t2(i,j).zeros();
        }
    }
    cout << "Successfully initialized t2 amplitudes." << endl;
}

double ccsolve::GetCoupledElement(int a, int b, int c, int d){
    double sm = 0.0;
    for(int i=0; i<nStates; i++){
        for(int j=0; j<nStates; j++){
            for(int k=0; k<nStates; k++){
                for(int l=0; l<nStates; l++){
                    sm += hfobject.C(a,i)*hfobject.C(b,j)*hfobject.C(c,k)*hfobject.C(d,l)*hfobject.coupledMatrix(i,j)(k,l);
                }
            }
        }
    }
    return sm;
}

double ccsolve::GetUncoupledElement(int a, int b){
    double sm = 0.0;
    for(int i=0; i<nStates; i++){
        for(int j=0; j<nStates; j++){
            sm += hfobject.C(a,i)*hfobject.C(b,j)*hfobject.Bs.h(i,j);
        }
    }
    return sm;
}

void ccsolve::SetupMinimizedBasis(){
    cout << "Setting up the minimized basis." << endl;
    vmin.set_size(nStates,nStates);
    for(int a=0; a<nStates; a++){
        for(int b=0; b<nStates; b++){
            vmin(a,b) = zeros<mat>(nStates,nStates);
            for(int c=0; c<nStates; c++){
                for(int d=0; d<nStates; d++){
                    vmin(a,b)(c,d) = GetCoupledElement(a,b,c,d);
                }
            }
        }
    }
    cout << "Finished setting up the minimized basis." << endl;
    //coupledMatrixMinimized.print();
}
