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
    //hfobject.C.print();
    cout << "CCSolve initialized." << endl;
    cout << "Identifying number of orbitals in the basis" << endl;
    cout << "Found " << nStates << " number of basis functions in the HF Rothaan expansion." << endl;
    //expandC();
    SetupMinimizedBasis();
    ExpandMinimizedBasis(); //Include spin orthogonality
    SetupT1();
    SetupT2();
    CCD();
    //fmin.print();
    cout << "Energy:" << energy() << endl;

    //hfobject.Bs.h.print();
}

void ccsolve::ExpandMinimizedBasis(){

    nStates*= 2;
    temp_mo = vmin;
    vmin.set_size(nStates, nStates);
    fmin.set_size(nStates, nStates);

    double val1 = 0.0;
    double val2 = 0.0;
    for(int a = 0; a<nStates; a++){
        for(int b = 0; b<nStates; b++){
            fmin(a,b) = GetUncoupledElement(a,b);
            vmin(a,b) = zeros(nStates,nStates);
            for(int i = 0; i<nStates; i++){
                for(int j=0; j<nStates; j++){

                    val1 = equalfunc(a%2,b%2) * equalfunc(i%2,j%2) * temp_mo(a/2,b/2)(i/2,j/2);
                    val2 = equalfunc(a%2,j%2) * equalfunc(i%2,b%2) * temp_mo(a/2,j/2)(i/2,b/2);
                    cout << val1 << " " << val2 << endl;
                    //cout << "Expanding the minimized basis." << endl;
                    vmin(a,b)(i,j) = val1 -val2;
                    //cout << "Expanding the minimized basis." << endl;
                }
            }
        }
    }
}

void ccsolve::expandC(){
    //expand the coefficient matrix to explicitly include spin dependence
    Cm.set_size(nStates,nStates);
    Cm.zeros();
    for(int i = 0; i<nStates/2; i++){
        for(int j = 0; j<nStates/2; j++){
            Cm(2*i,j) = hfobject.C(i,j);
            Cm(2*i+1,nStates/2+j) = hfobject.C(i,j);
        }
    }
    Cm.print();
}

void ccsolve::initT2(){
    for(int a=nElectrons; a<nStates; a++){
        for(int b=a+1; b<nStates; b++){
            for(int i=0;i<nElectrons;i++){
                for(int j=i+1;j<nElectrons;j++){
                    t2new(a,b)(i,j) = vmin(a,b)(i,j)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                }
            }
        }
    }
}

double ccsolve::CCDQ(int a, int b, int i, int j){
    double Qa = 0;
    double Qb = 0;
    double Qc = 0;
    double Qd = 0;
    //int a,b,i,j; //These should be included in the function call
    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qa += vmin(k,l)(c,d)*t2(c,d)(i,j)*t2(a,b)(k,l);
                }
            }
        }
    }
    //Qa/=4.0;

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qb += vmin(k,l)(c,d)*(t2(a,c)(i,k)*t2(b,d)(j,l) - t2(a,c)(j,k)*t2(b,d)(i,l));
                }
            }
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qc += -vmin(k,l)(c,d)*t2(d,c)(i,k)*t2(a,b)(l,j) + vmin(k,l)(c,d)*t2(d,c)(j,k)*t2(a,b)(l,i);
                }
            }
        }
    }
    //Qc *= .5;

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qd += -vmin(k,l)(c,d)*t2(a,c)(l,k)*t2(d,b)(i,j) + vmin(k,l)(c,d)*t2(b,c)(l,k)*t2(d,a)(i,j);
                }
            }
        }
    }
    //Qd *= .5;
    return 0.25*Qa + Qb + 0.5*(Qc + Qd);
}


double ccsolve::CCDL(int a, int b, int i, int j){
    double Da = 0.0;
    double Db = 0.0;
    double Dc = 0.0;
    double Dd = 0.0;
    double De = 0.0;
    double Df = 0.0;


    for(int c = nElectrons; c<nStates; c++){
        Da += fmin(b,c)*t2(a,c)(i,j) - fmin(a,c)*t2(b,c)(i,j);
    }

    for(int k = 0; k<nElectrons; k++){
        Db += -fmin(k,j)*t2(a,b)(i,k) + fmin(k,i)*t2(a,b)(j,k);
    }

    for(int c = nElectrons; c<nStates; c++){
        for(int d = nElectrons; d<nStates; d++){
            Dc += vmin(a,b)(c,d)*t2(c,d)(i,j);
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; k<nElectrons; k++){
            Dd += vmin(k,l)(i,j)*t2(a,b)(k,l);
        }
    }

    for(int k= 0; k<nElectrons; k++){
        for(int c = nElectrons; c<nStates; c++){
            De += vmin(k,b)(c,j)*t2(a,c)(i,k) - vmin(k,a)(c,j)*t2(b,c)(i,k) -vmin(k,b)(c,i)*t2(a,c)(j,k) + vmin(k,a)(c,i)*t2(b,c)(j,k);
        }
    }
    return 0*Da + 0*Db + .5*(Dc + Dd)  + De;

}

void ccsolve::CCD(){
    //Set up and solve the CCD-equation
    cout << "Performing CCD calculation." << endl;
    initT2(); //Set up initial guess following S-B, p.289

    cout << "Entering iterative scheme." << endl;
    //t2new.print();
    while(unconverged(.00001)){
        cout <<"Current energy: " << energy() << endl;
        t2 = t2new;
        //t2.print();
        double outfactored = 0.0;
        for(int a = nElectrons; a<nStates; a++){
            for(int b = a+1; b<nStates; b++){
                for(int i = 0; i<nElectrons; i++){
                    for(int j=i+1; j<nElectrons; j++){
                        outfactored = (-fmin(i,i) -fmin(j,j) + fmin(a,a) + fmin(b,b))*t2(a,b)(i,j);
                        //cout << "Denominator " << (-fmin(i,i) -fmin(j,j) + fmin(a,a) + fmin(b,b)) << i << j << a <<b << endl;
                        //cout << "Numerator "  << vmin(a,b)(i,j) + CCDL(a,b,i,j) + CCDQ(a,b,i,j)<<endl;
                        t2new(a,b)(i,j) = (vmin(a,b)(i,j) + CCDL(a,b,i,j) + CCDQ(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //What about the terms factored outside?
                    }
                }
            }
        }
    }
}

double ccsolve::energy(){
    double dE = 0;
    for(int i = 0; i<nElectrons; i++){
        for(int j = i+1; j<nElectrons; j++){
            for(int a=nElectrons; a<nStates; a++){
                for(int b=a+1; b<nStates; b++){
                    dE += vmin(i,j)(a,b)*t2new(a,b)(i,j);
                }
            }
        }
    }
    return dE/4.0;
}

double ccsolve::equalfunc(int a, int b){
    double ret = 0.0;
    if(a==b){
        ret = 1.0;
    }
    return ret;
}

bool ccsolve::unconverged(double tolerance){
    double diff = 0.0;
    bool condition = true;
    for(int p = 0; p<nStates; p++){
        for(int q = 0; q<nStates; q++){
            for(int r = 0; r<nStates; r++){
                for(int s = 0; s<nStates; s++){
                      diff += abs(t2(p,q)(r,s) - t2new(p,q)(r,s));
                }
            }
        }
    }
    if(diff<tolerance){condition = false;}
    if(diff!=diff){condition = false;}
    cout << "Difference in t2: " <<  diff << endl;
    return condition;
}

void ccsolve::SetupT1(){
    t1.set_size(nStates, nStates);
    t1.zeros();
    t1new.set_size(nStates, nStates);
    t1new.zeros();

    cout << "Successfully initialized t1 amplitudes." << endl;
}

void ccsolve::SetupT2(){
    t2.set_size(nStates, nStates);
    t2new.set_size(nStates, nStates);
    for(int i = 0; i<nStates; i++){
        for(int j =0;j<nStates; j++){
            t2(i,j).set_size(nStates, nStates);
            t2(i,j).zeros();
            t2new(i,j).set_size(nStates, nStates);
            t2new(i,j).zeros();
        }
    }
    cout << "Successfully initialized t2 amplitudes." << endl;
}

double ccsolve::GetCoupledElement(int a, int b, int c, int d){
    //this needs to be
    double sm = 0.0;
    for(int i=0; i<nStates; i++){
        for(int j=0; j<nStates; j++){
            for(int k=0; k<nStates; k++){
                for(int l=0; l<nStates; l++){

                    //sm += hfobject.C(a,i)*hfobject.C(b,j)*hfobject.C(c,k)*hfobject.C(d,l)*hfobject.coupledMatrix(i,j)(k,l);
                    sm += hfobject.C(a,i)*hfobject.C(b,j)*hfobject.C(c,k)*hfobject.C(d,l)*hfobject.coupledMatrix(i/2,j/2)(k/2,l/2);

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
    sm = 0.0;
    if(a==b){
        sm = hfobject.epsilon(a/2);
    }
    //return hfobject.Bs.h(a,b);
    return sm;
}

void ccsolve::SetupMinimizedBasis(){
    cout << "Setting up the minimized basis." << endl;
    vmin.set_size(nStates,nStates);
    fmin.set_size(nStates,nStates);
    for(int a=0; a<nStates; a++){
        for(int b=0; b<nStates; b++){

            fmin(a,b) = GetUncoupledElement(a,b);
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
