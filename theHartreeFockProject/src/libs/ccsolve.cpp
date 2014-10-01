#include "ccsolve.h"
#include <hartreefocksolver.h>

ccsolve::ccsolve()
{
}

ccsolve::ccsolve(hartreefocksolver object, int nElect)
{
    /* --- CCSolve Main ---
     * This class should be initialized with a hartreefocksolver object containing an energy-minimized basis for the Coupled-Cluster (CC) solver.
     * The CC Solver will then perform the following operations
     * 1. Set up a new basis from atomic to molecular orbitals
     * 2. Initialize the cluster amplitudes with a given initial value (Begin with CCD only)
     * 3. Solve the amplitude equations
     * 4. Return the CC-energy
    */
    hfobject = object;
    nElectrons = nElect; //Fermi level is defined by number of electrons
    nStates = hfobject.C.n_cols;
    cout << "CCSolve Found " << nStates << " number of basis functions in the RHF Roothaan expansion." << endl;
    SetupMinimizedBasis();
    ExpandMinimizedBasis(); //Include spin orthogonality
    SetupT1();
    SetupT2();
    CCSD();
    cout << "Energy:" << energy() << endl;
}

double ccsolve::GetUncoupledElement(int a, int b){
    double sm = 0;
    if(a==b){
        sm = hfobject.epsilon(a/2);
    }
    return sm;
}

double ccsolve::GetCoupledElement(int a, int b, int c, int d){
    //Unsure about this transformation
    double sm = 0.0;
    for(int i=0; i<nStates; i++){
        for(int j=0; j<nStates; j++){
            for(int k=0; k<nStates; k++){
                for(int l=0; l<nStates; l++){
                    //sm += hfobject.C(a,i)*hfobject.C(b,j)*hfobject.C(k,c)*hfobject.C(l,d)*hfobject.coupledMatrix(i/2,j/2)(k/2,l/2);
                    //sm += hfobject.C(a,i)*hfobject.C(b,j)*hfobject.C(c,k)*hfobject.C(d,l)*hfobject.Bs.v(i,j)(k,l);
                    sm += hfobject.C(i,a)*hfobject.C(j,b)*hfobject.C(k,c)*hfobject.C(l,d)*hfobject.Bs.v(i,j)(k,l);
                }
            }
        }
    }
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
                    //vmin(a,b)(c,d) = hfobject.P(b,c)*hfobject.P(a,d)*hfobject.Bs.v(a,b)(c,d);
                    vmin(a,b)(c,d) = GetCoupledElement(a,b,c,d);
                }
            }
        }
    }
    cout << "Finished setting up the minimized basis." << endl;
}

void ccsolve::ExpandMinimizedBasis(){
    nStates*= 2;
    temp_mo = vmin;
    vmin.set_size(nStates, nStates);
    fmin.set_size(nStates, nStates);

    double val1 = 0.0;
    double val2 = 0.0;
    for(int p = 0; p<nStates; p++){
        for(int q = 0; q<nStates; q++){
            //previously aibj
            fmin(p,q) = GetUncoupledElement(p,q);
            vmin(p,q) = zeros(nStates,nStates);
            for(int r = 0; r<nStates; r++){
                for(int s=0; s<nStates; s++){

                    val1 = equalfunc(p%2,q%2) * equalfunc(r%2,s%2) * temp_mo(p/2,q/2)(r/2,s/2);
                    //val2 = equalfunc(a%2,j%2) * equalfunc(i%2,b%2) * temp_mo(a/2,j/2)(i/2,b/2); //Originals

                    val2 = equalfunc(p%2,r%2) * equalfunc(q%2,s%2) * temp_mo(p/2,q/2)(s/2,r/2); //Originals

                    //vmin(a,i)(b,j) = val1 -val2; //original

                    //vmin(a,b)(i,j) = val1 -val2;
                    vmin(p,q)(r,s) = hfobject.Bs.state(p,q,r,s, temp_mo(p/2,q/2)(r/2,s/2), temp_mo(p/2,q/2)(s/2,r/2));
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

double ccsolve::CCSD_Single(int a,int i){
    double S1 = 0;
    double S2a = 0;
    double S2b = 0;
    double S2c = 0;
    double S3a = 0;
    double S3b = 0;
    double S3c = 0;
    double S4a = 0;
    double S4b = 0;
    double S4c = 0;
    double S5a = 0;
    double S5b = 0;
    double S5c = 0;
    double S6 = 0;
    int k, l, c, d;

    S1 += fmin(a,i);

    for(k=0;k<nElectrons;k++){
        for(c=nElectrons;c<nStates;c++){
            S2a+=fmin(k,c)*t2(a,c)(i,k);
            for(d=nElectrons;d<nStates;d++){
                S2b+=vmin(a,k)(c,d)*t2(c,d)(i,k);
            }
        }
    }


    for(k=0;k<nElectrons;k++){
        for(l=0;l<nElectrons;l++){
            for(c=nElectrons;c<nStates;c++){
                S2c-=vmin(k,l)(i,c)*t2(a,c)(k,l);
            }
        }
    }

    for(c=nElectrons;c<nStates;c++){
        S3a+=fmin(a,c)*t1(c,i);
    }


    for(k=0;k<nElectrons;k++){
        S3b-=fmin(k,i)*t1(a,k);
        for(c=nElectrons;c<nStates;c++){
            S3c += vmin(a,k)(i,c)*t1(c,k);
        }
    }

    for(k=0;k<nElectrons;k++){
        for(l=0;l<nElectrons;l++){
            for(c=nElectrons;c<nStates;c++){
                S5c-=vmin(k,l)(i,c)*t1(a,k)*t1(c,l);
                for(d=nElectrons;d<nStates;d++){
                    S4a-=vmin(k,l)(c,d)*t1(c,i)*t2(a,d)(k,l);
                    S4b-=vmin(k,l)(c,d)*t1(a,k)*t2(c,d)(i,l);
                    S4c+=vmin(k,l)(c,d)*t1(c,k)*t2(d,a)(l,i);
                    S6-=vmin(k,l)(c,d)*t1(c,i)*t1(a,k)*t1(d,l);
                }
            }
        }
    }

    for(k=0;k<nElectrons;k++){
        for(c=nElectrons;c<nStates;c++){
            S5a -= fmin(k,c)*t1(c,i)*t1(a,k);
            for(d=nElectrons; d<nStates; d++){
                S5b += vmin(a,k)(c,d)*t1(c,i)*t1(d,k);
            }
        }
    }

    return S1+S2a+.5*S2b+.5*S2c+   0*S3a+0*S3b  +S3c+.5*S4a+.5*S4b+S4c+S5a+S5b+S5c+S6;
}

double ccsolve::CCSD_Double(int a, int b, int i, int j){
    int k,l,c,d;
    int nE = nElectrons;
    int nS = nStates;
    double D4a = 0;
    double D4b = 0;

    double D5a = 0;
    double D5b = 0;
    double D5c = 0;
    double D5d = 0;
    double D5e = 0;
    double D5f = 0;
    double D5g = 0;
    double D5h = 0;

    double D6a = 0;
    double D6b = 0;
    double D6c = 0;

    double D7a = 0;
    double D7b = 0;
    double D7c = 0;
    double D7d = 0;
    double D7e = 0;

    double D8a = 0;
    double D8b = 0;

    double D9  = 0;


    for(c=nE;c<nS;c++){
        D4a += vmin(a,b)(c,j)*t1(c,i) - vmin(a,b)(c,i)*t1(c,j);
    }

    for(k=0;k<nE;k++){
        D4b += -vmin(k,b)(i,j)*t1(a,k) + vmin(k,a)(i,j)*t1(b,k);
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            D5a += -fmin(k,c)*t1(c,i)*t2(a,b)(k,j) + fmin(k,c)*t1(c,j)*t2(a,b)(k,i);
            D5b += -fmin(k,c)*t1(a,k)*t2(c,b)(i,j) + fmin(k,c)*t1(b,k)*t2(c,a)(i,j);
            for(d=nE;d<nS;d++){
                D5c += vmin(a,k)(c,d)*t1(c,i)*t2(d,b)(k,j);
                D5c -= vmin(a,k)(c,d)*t1(c,j)*t2(d,b)(k,i);
                D5c -= vmin(b,k)(c,d)*t1(c,i)*t2(d,a)(k,j);
                D5c += vmin(b,k)(c,d)*t1(c,j)*t2(d,a)(k,i);
            }
        }
    }

    for(k=0;k<nE;k++){
        for(l=0;l<nE;l++){
            for(c=nE;c<nS;c++){
                D5d -= vmin(k,l)(i,c)*t1(a,k)*t2(c,b)(l,j);
                D5d += vmin(k,l)(j,c)*t1(a,k)*t2(c,b)(l,i);
                D5d += vmin(k,l)(i,c)*t1(b,k)*t2(c,a)(l,j);
                D5d -= vmin(k,l)(j,c)*t1(b,k)*t2(c,a)(l,i);
            }
        }
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            for(d=nE;d<nS;d++){
                D5e -= vmin(k,b)(c,d)*t1(a,k)*t2(c,d)(i,j);
                D5e += vmin(k,a)(c,d)*t1(b,k)*t2(c,d)(i,j);
            }
        }
    }

    for(k=0;k<nE;k++){
        for(l=0;l<nE;l++){
            for(c=nE;c<nS;c++){
                D5f += vmin(k,l)(c,j)*t1(c,i)*t2(a,b)(k,l);
                D5f -= vmin(k,l)(c,i)*t1(c,j)*t2(a,b)(k,l);
            }
        }
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            for(d=nE;d<nS;d++){
                D5g += vmin(k,a)(c,d)*t1(c,k)*t2(d,b)(i,j);
                D5g -= vmin(k,b)(c,d)*t1(c,k)*t2(d,a)(i,j);
            }
        }
    }


    for(k=0;k<nE;k++){
        for(l=0;l<nE;l++){
            for(c=nE;c<nS;c++){
                D5h -= vmin(k,l)(c,i)*t1(c,k)*t2(a,b)(l,j);
                D5h += vmin(k,l)(c,j)*t1(c,k)*t2(a,b)(l,i);
            }
        }
    }

    for(c=nE;c<nS;c++){
        for(d=nE;d<nS;d++){
            D6a += vmin(a,b)(c,d)*t1(c,i)*t1(d,j);
        }
    }

    for(k=0;k<nE;k++){
        for(l=0;l<nE;l++){
            D6b += vmin(k,l)(i,j)*t1(a,k)*t1(b,l);
        }
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            D6c -= vmin(k,b)(c,j)*t1(c,i)*t1(a,k);
            D6c += vmin(k,b)(c,i)*t1(c,j)*t1(a,k);
            D6c += vmin(k,a)(c,j)*t1(c,i)*t1(b,k);
            D6c -= vmin(k,a)(c,i)*t1(c,j)*t1(b,k);
        }
    }

    for(k=0;k<nE;k++){
        for(l=0;l<nE;l++){
            for(c=nE;c<nS;c++){
                for(d=nE; d<nS;d++){
                    D7a += vmin(k,l)(c,d)*t1(c,i)*t1(d,j)*t2(a,b)(k,l);
                    D7b += vmin(k,l)(c,d)*t1(a,k)*t1(b,l)*t2(c,d)(i,j);

                    D7c -= vmin(k,l)(c,d)*t1(c,i)*t1(a,k)*t2(d,b)(l,j);
                    D7c += vmin(k,l)(c,d)*t1(c,j)*t1(a,k)*t2(d,b)(l,i);
                    D7c += vmin(k,l)(c,d)*t1(c,i)*t1(b,k)*t2(d,a)(l,j);
                    D7c -= vmin(k,l)(c,d)*t1(c,j)*t1(b,k)*t2(d,a)(l,i);

                    D7d -= vmin(k,l)(c,d)*t1(c,k)*t1(d,i)*t2(a,b)(l,j);
                    D7d += vmin(k,l)(c,d)*t1(c,k)*t1(d,j)*t2(a,b)(l,i);

                    D7e -= vmin(k,l)(c,d)*t1(c,k)*t1(a,l)*t2(d,b)(i,j);
                    D7e += vmin(k,l)(c,d)*t1(c,k)*t1(b,l)*t2(d,a)(i,j);
                }
            }
        }
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            for(d=nE;d<nS;d++){
                D8a += vmin(k,b)(c,d)*t1(c,i)*t1(a,k)*t1(d,j);
                D8a -= vmin(k,a)(c,d)*t1(c,i)*t1(b,k)*t1(d,j);
            }
        }
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            for(d=nE;d<nS;d++){
                D8b += vmin(k,l)(c,j)*t1(c,i)*t1(a,k)*t1(b,l);
                D8b -= vmin(k,l)(c,i)*t1(c,j)*t1(a,k)*t1(b,l);
            }
        }
    }

    for(k=0;k<nE;k++){
        for(l=0;l<nE;l++){
            for(c=nE;c<nS;c++){
                for(d=nE; d<nS;d++){
                    D9 += vmin(k,l)(c,d)*t1(c,i)*t1(d,j)*t1(a,k)*t1(b,l);
                }
            }
        }
    }

    return D4a+D4b+D5a+D5b+D5c+.5*D5e+D5g+.5*D5f+D5h + D6a + D6b + D6c + .5*D7a + .5*D7b + D7c + D7d + D7e + D8a + D8b + D9;
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
                    Qb += vmin(k,l)(c,d)*(t2(a,c)(i,k)*t2(d,b)(j,l) - t2(a,c)(j,k)*t2(d,b)(i,l));
                }
            }
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qc -= vmin(k,l)(c,d)*t2(d,c)(i,k)*t2(a,b)(l,j);
                    Qc += vmin(k,l)(c,d)*t2(d,c)(j,k)*t2(a,b)(l,j);
                }
            }
        }
    }
    //Qc *= .5;

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qd -= vmin(k,l)(c,d)*t2(a,c)(l,k)*t2(d,b)(i,j);
                    Qd += vmin(k,l)(c,d)*t2(b,c)(l,k)*t2(d,a)(i,j);
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
        Db -= fmin(k,j)*t2(a,b)(i,k);
        Db += fmin(k,i)*t2(a,b)(j,k);
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
            De += vmin(k,b)(c,j)*t2(a,c)(i,k);
            De -= vmin(k,a)(c,j)*t2(b,c)(i,k);
            De -= vmin(k,b)(c,i)*t2(a,c)(j,k);
            De += vmin(k,a)(c,i)*t2(b,c)(j,k);
        }
    }
    return 0*Da + 0*Db + .5*(Dc + Dd)  + De;

}

void ccsolve::CCSD(){
    eprev = 0.0;
    //Set up and solve the CCSD-equation
    cout << "Performing CCSD calculation." << endl;
    initT1(); //initial guess for the amplitudes
    initT2(); //Set up initial guess following S-B, p.289
    while(unconverged(.00000001)){
        cout <<"Current energy: " << energy() << endl;
        t2 = t2new; //updating the amplitudes
        t1 = t1new;
        eprev = energy(); //updating the energy
        double outfactored = 0.0;
        for(int b = nElectrons; b<nStates; b++){
            for(int a = b+1; a<nStates; a++){
                for(int j = 0; j<nElectrons; j++){
                    for(int i=j+1; i<nElectrons; i++){
                        //outfactored = (-fmin(i,i) -fmin(j,j) + fmin(a,a) + fmin(b,b))*t2(a,b)(i,j);
                        //t2new(b,a)(j,i) = (vmin(a,b)(i,j) + CCSD_Double(i,j,a,b) +  CCDL(a,b,i,j) + CCDQ(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //What about the terms factored outside?
                        t2new(a,b)(i,j) = (vmin(a,b)(i,j) + CCSD_Double(a,b,i,j) +  CCDL(a,b,i,j) + CCDQ(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));

                    }
                }
            }
        }
        for(int a=nElectrons;a<nStates;a++){
            for(int i=0;i<nElectrons;i++){
                t1new(a,i) = CCSD_Single(a,i)/(fmin(i,i)-fmin(a,a));

            }
        }
        cout << t1new.max() << endl;
    }
}

void ccsolve::CCD(){
    eprev = 0.0;
    //Set up and solve the CCSD-equation
    cout << "Performing CCD calculation." << endl;
    //initT1(); //initial guess for the amplitudes
    initT2(); //Set up initial guess following S-B, p.289
    while(unconverged(.00000001)){
        cout <<"Current energy: " << energy() << endl;
        t2 = t2new; //updating the amplitudes
        //t1 = t1new;
        eprev = CCDenergy(); //updating the energy
        for(int b = nElectrons; b<nStates; b++){
            for(int a = b+1; a<nStates; a++){
                for(int j = 0; j<nElectrons; j++){
                    for(int i=j+1; i<nElectrons; i++){
                        t2new(b,a)(j,i) = (vmin(a,b)(i,j) + CCDL(a,b,i,j) + CCDQ(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //What about the terms factored outside?
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
                    dE += .25*vmin(i,j)(a,b)*t2new(a,b)(i,j);
                    dE += .5*vmin(i,j)(a,b)*t1new(a,i)*t1new(b,j);
                }
            }
        }
    }

    for(int i = 0; i<nElectrons; i++){
        for(int a=nElectrons; a<nStates; a++){
            dE += fmin(i,a)*t1new(a,i);

        }
    }

    return dE;
}


double ccsolve::equalfunc(int a, int b){
    double ret = 0.0;
    if(a==b){
        ret = 1.0;
    }
    return ret;
}


double ccsolve::CCDenergy(){
    double dE = 0;
    for(int i = 0; i<nElectrons; i++){
        for(int j = i+1; j<nElectrons; j++){
            for(int a=nElectrons; a<nStates; a++){
                for(int b=a+1; b<nStates; b++){
                    dE += .25*vmin(i,j)(a,b)*t2new(a,b)(i,j);
                }
            }
        }
    }
    return dE;
}

bool ccsolve::unconverged(double tolerance){
    bool condition = true;
    if(abs(energy()-eprev)<tolerance){
        cout << "Converged at " << abs(energy()-eprev) << endl;
        condition = false;
    }
    if(energy()!=energy()){condition = false;}
    return condition;
}

void ccsolve::initT2(){
    for(int a=nElectrons; a<nStates; a++){
        for(int b=nElectrons; b<nStates; b++){
            for(int i=0;i<nElectrons;i++){
                for(int j=0;j<nElectrons;j++){
                    t2new(a,b)(i,j) = vmin(a,b)(i,j)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                }
            }
        }
    }
}

void ccsolve::initT1(){
    for(int a=nElectrons; a<nStates; a++){
        for(int i=0; i<nElectrons; i++){
            t1new(a,i) = fmin(a,i)/(fmin(i,i) - fmin(a,a));
        }
    }
}

void ccsolve::SetupT1(){
    t1.set_size(nStates, nStates);
    t1.zeros();
    t1new.set_size(nStates, nStates);
    t1new.zeros();

    cout << "Initialized t1 amplitudes." << endl;
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
    cout << "Initialized t2 amplitudes." << endl;
}



