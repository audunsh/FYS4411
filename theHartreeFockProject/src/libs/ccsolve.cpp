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

    //hfobject.C.col(3)*=-1;
    //hfobject.C.col(5)*=-1;

    cout << "CCSolve Found " << nStates << " number of basis functions in the RHF Roothaan expansion." << endl;
    SetupMinimizedBasis();
    ExpandMinimizedBasis(); //Include spin orthogonality
    SetupT1();
    SetupT2();
    //retranslate();
    CCD();
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
                    //sm += hfobject.C(a,i)*hfobject.C(b,j)*hfobject.C(c,k)*hfobject.C(d,l)*hfobject.coupledMatrix(i,j)(k,l);
                    //sm += hfobject.C(a,i)*hfobject.C(b,j)*hfobject.C(c,k)*hfobject.C(d,l)*hfobject.Bs.v(i,j)(k,l);

                    sm += hfobject.C(i,a)*hfobject.C(j,b)*hfobject.C(k,c)*hfobject.C(l,d)*hfobject.Bs.v(i,j)(k,l);
                    //sm += hfobject.C(i,a)*hfobject.C(j,b)*hfobject.C(k,c)*hfobject.C(l,d)*hfobject.coupledMatrix(i,k)(j,l);
                    //sm += hfobject.C(i,a)*hfobject.C(j,b)*hfobject.C(k,c)*hfobject.C(l,d)*hfobject.coupledMatrix(i,j)(k,l);
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
    fmin.zeros();
    for(int a=0; a<nStates; a++){
        for(int b=0; b<nStates; b++){
            //fmin(a,b) = GetUncoupledElement(a,b);
            vmin(a,b) = zeros<mat>(nStates,nStates);}}
    for(int a=0; a<nStates; a++){
        //fmin(a,a) = hfobject.epsilon(a/2);
        for(int b=0; b<nStates; b++){
            for(int c=0; c<nStates; c++){
                for(int d=0; d<nStates; d++){
                    //vmin(a,b)(c,d) = hfobject.P(b,c)*hfobject.P(a,d)*hfobject.Bs.v(a,b)(c,d);
                    vmin(a,b)(c,d) = GetCoupledElement(a,b,c,d) ; //Adjusting this parameter to compensate for notational differences
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
    fmin.zeros();
    //hfobject.C.col(4)*=-1;

    double val1 = 0.0;
    double val2 = 0.0;
    for(int p = 0; p<nStates; p++){
        for(int q = 0; q<nStates; q++){
            //previously aibj
            fmin(p,q) = GetUncoupledElement(p,q);
            vmin(p,q) = zeros(nStates,nStates);
        }
    }
    for(int p = 0; p<nStates; p++){
        for(int q = 0; q<nStates; q++){
            for(int r = 0; r<nStates; r++){
                for(int s=0; s<nStates; s++){

                    val1 = equalfunc(p%2,q%2) * equalfunc(r%2,s%2) * temp_mo(p/2,q/2)(r/2,s/2);
                    val2 = equalfunc(p%2,s%2) * equalfunc(r%2,q%2) * temp_mo(p/2,s/2)(r/2,q/2);
                    vmin(p,r)(q,s) = val1 - val2;



                    //vmin(p,r)(q,s) = hfobject.Bs.state(p,r,q,s, temp_mo(p/2,q/2)(r/2,s/2), temp_mo(p/2,s/2)(r/2,q/2)); //THIS PRODUCES SAME RESULTS (vmin) AS FROM NORDLI
                    //vmin(p,r)(q,s) = hfobject.Bs.state(p,r,q,s, temp_mo(p/2,r/2)(q/2,s/2), temp_mo(p/2,r/2)(q/2,s/2)); //THIS PRODUCES SAME RESULTS (vmin) AS FROM NORDLI

                    //vmin(p,q)(r,s) = hfobject.Bs.state(p,q,r,s, temp_mo(p/2,q/2)(r/2,s/2), temp_mo(p/2,q/2)(s/2,r/2)); //Trying to make minor changes in indexing
                    //vmin(p,q)(r,s) = hfobject.Bs.state(p,q,r,s, temp_mo(p/2,q/2)(r/2,s/2), temp_mo(p/2,q/2)(s/2,r/2)); //Minor editing, producing errors
                }
            }
        }
    }
    cout << endl;
    //vmin.print();
    // vmin(13,12);
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

    //return S1+S2a+.5*S2b+.5*S2c+   0*S3a+0*S3b  +S3c+.5*S4a+.5*S4b+S4c+S5a+S5b+S5c+S6;
    return .5*S2b+.5*S2c+   0*S3a+0*S3b  +S3c+.5*S4a+.5*S4b+S4c+S5b+S5c+S6; //HF-case
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
        D4a += vmin(a,b)(c,j)*t1(c,i);
        D4a -= vmin(a,b)(c,i)*t1(c,j);
    }

    for(k=0;k<nE;k++){
        D4b -= vmin(k,b)(i,j)*t1(a,k);
        D4b += vmin(k,a)(i,j)*t1(b,k);
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            D5a -= fmin(k,c)*t1(c,i)*t2(a,b)(k,j);
            D5a += fmin(k,c)*t1(c,j)*t2(a,b)(k,i);

            D5b -= fmin(k,c)*t1(a,k)*t2(c,b)(i,j);
            D5b += fmin(k,c)*t1(b,k)*t2(c,a)(i,j);

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
    //Quadratic contributions to the CCD amplitude equation
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
                    //Qb += vmin(k,l)(c,d)*t2(a,c)(i,k)*t2(d,b)(l,j); (*) Conflictin notation in SB
                    //Qb -= vmin(k,l)(c,d)*t2(a,c)(j,k)*t2(d,b)(l,i);

                    Qb += vmin(k,l)(c,d)*t2(a,c)(i,k)*t2(b,d)(j,l);
                    Qb -= vmin(k,l)(c,d)*t2(a,c)(j,k)*t2(b,d)(i,l);
                }
            }
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qc -= vmin(k,l)(c,d)*t2(c,d)(k,i)*t2(a,b)(l,j);
                    Qc += vmin(k,l)(c,d)*t2(c,d)(k,j)*t2(a,b)(l,i);
                }
            }
        }
    }
    //Qc *= .5;

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qd -= vmin(k,l)(c,d)*t2(c,a)(k,l)*t2(d,b)(i,j);
                    Qd += vmin(k,l)(c,d)*t2(c,b)(k,l)*t2(d,a)(i,j);
                }
            }
        }
    }
    //Qd *= .5;
    return 0.25*Qa + Qb + 0.5*(Qc + Qd);
}

double ccsolve::CCDQ2(int a, int b, int i, int j){
    //Quadratic contributions to the CCD amplitude equation
    //Call the elements in the following order: vmin(i,j)(a,b)
    double Qa = 0;
    double Qb = 0;
    double Qc = 0;
    double Qd = 0;
    //int a,b,i,j; //These should be included in the function call
    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qa += vmin(k,l)(c,d)*t2p(c,d)(i,j)*t2p(a,b)(k,l);
                }
            }
        }
    }
    //Qa/=4.0;

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    //Qb += vmin(k,l)(c,d)*t2(a,c)(i,k)*t2(d,b)(l,j); (*) Conflictin notation in SB
                    //Qb -= vmin(k,l)(c,d)*t2(a,c)(j,k)*t2(d,b)(l,i);

                    Qb += vmin(k,l)(c,d)*t2p(a,c)(i,k)*t2p(b,d)(j,l);
                    Qb -= vmin(k,l)(c,d)*t2p(a,c)(j,k)*t2p(b,d)(i,l);
                }
            }
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qc -= vmin(k,l)(c,d)*t2p(d,c)(i,k)*t2p(a,b)(l,j);
                    Qc += vmin(k,l)(c,d)*t2p(d,c)(j,k)*t2p(a,b)(l,i);
                }
            }
        }
    }
    //Qc *= .5;

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qd -= vmin(k,l)(c,d)*t2p(a,c)(l,k)*t2p(d,b)(i,j);
                    Qd += vmin(k,l)(c,d)*t2p(b,c)(l,k)*t2p(d,a)(i,j);
                }
            }
        }
    }
    //Qd *= .5;
    return 0.25*Qa + Qb + 0.5*(Qc + Qd);
}

double ccsolve::CCDL2(int a, int b, int i, int j){
    //Linear contributions to the CCD amplitude equation.
    //Call the elements in the following order: v(i,j)(a,b)
    double L1a = 0.0;
    double L1b = 0.0;
    double L2a = 0.0;
    double L2b = 0.0;
    double L2c = 0.0;


    for(int c = nElectrons; c<nStates; c++){
        L1a += fmin(b,c)*t2c(a,c)(i,j) - fmin(a,c)*t2c(b,c)(i,j);
    }

    for(int k = 0; k<nElectrons; k++){
        L1b -= fmin(k,j)*t2c(a,b)(i,k);
        L1b += fmin(k,i)*t2c(a,b)(j,k);
    }

    for(int c = nElectrons; c<nStates; c++){
        for(int d = nElectrons; d<nStates; d++){
            L2a += vmin(a,b)(c,d)*t2c(c,d)(i,j); //H2O is highly sensitive to this contribution (2/10 2014)

            //L2a += vmin(a,c)(b,d)*t2(c,d)(i,j); //Shuffling some indices
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; k<nElectrons; k++){
            //L2b += vmin(i,j)(k,l)*t2(a,b)(k,l); //(*) Conflicting notation in SB
            L2b += vmin(i,j)(k,l)*t2c(a,b)(k,l);

            //L2b += vmin(k,i)(l,j)*t2(a,b)(k,l); //shuffling some indices
        }
    }

    for(int k= 0; k<nElectrons; k++){
        for(int c = nElectrons; c<nStates; c++){
            L2c += vmin(k,b)(c,j)*t2c(a,c)(i,k);
            L2c -= vmin(k,a)(c,j)*t2c(b,c)(i,k);
            L2c -= vmin(k,b)(c,i)*t2c(a,c)(j,k);
            L2c += vmin(k,a)(c,i)*t2c(b,c)(j,k);



            //L2c += vmin(b,k)(c,j)*t2c(a,c)(i,k);
            //L2c -= vmin(b,k)(c,i)*t2c(a,c)(j,k);
            //L2c -= vmin(a,k)(c,j)*t2c(b,c)(i,k);
            //L2c += vmin(a,k)(c,i)*t2c(b,c)(j,k);
        }
    }
    return .5*(L2a + L2b) + L2c; //0*.65*L2c;

}

double ccsolve::CCDL(int a, int b, int i, int j){
    //Linear contributions to the CCD amplitude equation.
    double L1a = 0.0;
    double L1b = 0.0;
    double L2a = 0.0;
    double L2b = 0.0;
    double L2c = 0.0;


    for(int c = nElectrons; c<nStates; c++){
        L1a += fmin(b,c)*t2(a,c)(i,j) - fmin(a,c)*t2(b,c)(i,j);
    }

    for(int k = 0; k<nElectrons; k++){
        L1b -= fmin(k,j)*t2(a,b)(i,k);
        L1b += fmin(k,i)*t2(a,b)(j,k);
    }

    for(int c = nElectrons; c<nStates; c++){
        for(int d = nElectrons; d<nStates; d++){
            L2a += vmin(a,b)(c,d)*t2(c,d)(i,j); //H2O is highly sensitive to this contribution (2/10 2014)

            //L2a += vmin(a,c)(b,d)*t2(c,d)(i,j); //Shuffling some indices
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; k<nElectrons; k++){
            //L2b += vmin(i,j)(k,l)*t2(a,b)(k,l); //(*) Conflicting notation in SB
            L2b += vmin(k,l)(i,j)*t2(a,b)(k,l);

            //L2b += vmin(k,i)(l,j)*t2(a,b)(k,l); //shuffling some indices
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int c = nElectrons; c<nStates; c++){
            L2c -= vmin(a,k)(c,j)*t2(c,b)(i,k);
            L2c += vmin(a,k)(c,i)*t2(c,b)(j,k);
            L2c += vmin(b,k)(c,j)*t2(c,a)(i,k);
            L2c -= vmin(b,k)(c,i)*t2(c,a)(j,k);

            //L2c += vmin(k,b)(c,j)*t2(a,c)(i,k);
            //L2c -= vmin(k,b)(c,i)*t2(a,c)(j,k);
            //L2c -= vmin(k,a)(c,j)*t2(b,c)(i,k);
            //L2c += vmin(k,a)(c,i)*t2(b,c)(j,k);

        }
    }
    //return 0*(L1a + L1b) + .5*(L2a + L2b)  + L2c; //0*.65*L2c;
    return 0*L1a + 0*L1b + .5*L2a + .5*L2b  + L2c; //0*.65*L2c;

}

void ccsolve::CCSD(){
    eprev = 0.0;
    maxiter = 100;
    counter = 0;
    //Set up and solve the CCSD-equation

    cout << "Performing CCSD calculation." << endl;
    initT1(); //initial guess for the amplitudes
    initT2(); //Set up initial guess following S-B, p.289
    t20 = t2new;

    while(unconverged(.00000001)){

        t2 = t2new; //updating the amplitudes
        t1 = t1new;
        cout <<"::Current energy: " << energy() << endl;
        eprev = energy(); //updating the energy

        for(int b = nElectrons; b<nStates; b++){
            for(int a = b+1; a<nStates; a++){
                for(int j = 0; j<nElectrons; j++){
                    for(int i=j+1; i<nElectrons; i++){
                        //t2new(a,b)(i,j) = (vmin(i,j)(a,b) + CCSD_Double(a,b,i,j) +  CCDL(a,b,i,j) + CCDQ(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //Working, original
                        t2new(a,b)(i,j) = (vmin(a,b)(i,j) + CCSD_Double(a,b,i,j) +  CCDL(a,b,i,j) + CCDQ(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //Faulty copy


                        t2new(b,a)(i,j) = -t2new(a,b)(i,j);
                        t2new(a,b)(j,i) = -t2new(a,b)(i,j);
                        t2new(b,a)(j,i) =  t2new(a,b)(i,j);
                    }
                }
            }
        }
        for(int a=nElectrons;a<nStates;a++){
            for(int i=0;i<nElectrons;i++){
                t1new(a,i) = CCSD_Single(a,i)/(fmin(i,i)-fmin(a,a));
            }
        }
    }
}

void ccsolve::retranslate(){
    temp_mo = vmin;
    int a,b,c,d;
    for(a=0;a<nStates;a++){
        for(b=0;b<nStates;b++){
            for(c=0;c<nStates;c++){
                for(d=0;d<nStates;d++){
                    vmin(a,b)(c,d) = temp_mo(a,d)(c,b);
                }
            }
        }
    }
}

void ccsolve::CCD(){
    //retranslate();
    eprev = 4.0;
    maxiter = 100;
    counter = 0;
    //Set up and solve the CCSD-equation
    cout << "Performing CCD calculation." << endl;
    //initT1(); //initial guess for the amplitudes
    initT2(); //Set up initial guess following S-B, p.289
    //t20 = t2new;
    while(unconverged(.00000001)){
        cout <<"Current energy: " << CCDenergy() << endl;


        //t2 = t2new; //updating the amplitudes
        //t1 = t1new;
        eprev = CCDenergy(); //updating the energy
        for(int b = nElectrons; b<nStates; b++){
            for(int a = b+1; a<nStates; a++){
                for(int j = 0; j<nElectrons; j++){
                    for(int i=j+1; i<nElectrons; i++){

                        //t2new(a,b)(i,j) = (vmin(i,j)(a,b) + CCDL(a,b,i,j) + CCDQ(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //Original, working until 1. iteration
                        //t2new(a,b)(i,j) = (vmin(i,j)(a,b) + CCDL2(a,b,i,j) + CCDQ2(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        //t2new(a,b)(i,j) = (vmin(a,b)(i,j) + CCDL2(a,b,i,j) + 0*CCDQ2(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        t2new(a,b)(i,j) = (vmin(a,b)(i,j) + CCDL2(a,b,i,j) + CCDQ2(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));

                        t2new(b,a)(i,j) = -t2new(a,b)(i,j);
                        t2new(a,b)(j,i) = -t2new(a,b)(i,j);
                        t2new(b,a)(j,i) =  t2new(a,b)(i,j);

                        //t2new(i,j)(a,b) = (vmin(a,b)(i,j) + CCDL(a,b,i,j) + CCDQ(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //Original, working until 1. iteration
                        //t2new(i,j)(b,a) = -t2new(i,j)(a,b);
                        //t2new(j,i)(a,b) = -t2new(i,j)(a,b);
                        //t2new(j,i)(b,a) =  t2new(i,j)(a,b);
                    }
                }
            }
        }

        t2p = t2c;
        t2c = t2new;
        //t2p = t2new;

    }
}


double ccsolve::energy(){
    double dE1 = 0;
    double dE2 = 0;
    double dE3 = 0;
    for(int i = 0; i<nElectrons; i++){
        for(int j = 0; j<nElectrons; j++){
            for(int a=nElectrons; a<nStates; a++){
                for(int b=nElectrons; b<nStates; b++){
                    //dE += .25*vmin(i,j)(a,b)*t2new(a,b)(i,j);
                    //dE +=  .5*vmin(i,j)(a,b)*t1new(a,i)*t1new(b,j); //originals, working

                    //dE += .25*vmin(i,j)(a,b)*t2new(a,b)(i,j);
                    //dE +=  .5*vmin(i,j)(a,b)*t1new(a,i)*t1new(b,j);

                    dE1 += vmin(i,j)(a,b)*t2new(a,b)(i,j);
                    dE2 += vmin(i,j)(a,b)*t1new(a,i)*t1new(b,j);
                }
            }
        }
    }

    for(int i = 0; i<nElectrons; i++){
        for(int a=nElectrons; a<nStates; a++){
            //since i!=a this will never contribute to the energy when using a HF SD.
            dE3 += fmin(i,a)*t1new(a,i);
        }
    }

    return .25*dE1 + .5*dE2 + dE3;
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
        for(int j = 0; j<nElectrons; j++){
            for(int a=nElectrons; a<nStates; a++){
                for(int b=nElectrons; b<nStates; b++){
                    //dE += vmin(i,j)(a,b)*t2new(a,b)(i,j);
                    dE += .25*vmin(i,j)(a,b)*t2c(a,b)(i,j);


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
    if(counter > maxiter){
        cout << "Maximum number of iterations exceeded. (" << maxiter << ")" << endl;
        condition = false;}
    counter += 1;
    return condition;
}


void ccsolve::initT2(){
    for(int a=nElectrons; a<nStates; a++){
        for(int b=nElectrons; b<nStates; b++){
            for(int i=0;i<nElectrons;i++){
                for(int j=0;j<nElectrons;j++){
                    //t2new(a,b)(i,j) = vmin(i,j)(a,b)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //Original, working
                    t2c(a,b)(i,j) = vmin(a,b)(i,j)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));

                }
            }
        }
    }
}


void ccsolve::initT1(){
    for(int a=nElectrons; a<nStates; a++){
        for(int i=0; i<nElectrons; i++){
            t1new(i,a) = fmin(a,i)/(fmin(i,i) - fmin(a,a));
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
    t2c.set_size(nStates, nStates);
    t2p.set_size(nStates, nStates);
    t2new.set_size(nStates, nStates);
    for(int i = 0; i<nStates; i++){
        for(int j =0;j<nStates; j++){
            t2(i,j).set_size(nStates, nStates);
            t2(i,j).zeros();
            t2new(i,j).set_size(nStates, nStates);
            t2new(i,j).zeros();
            t2c(i,j).set_size(nStates, nStates);
            t2c(i,j).zeros();
            t2p(i,j).set_size(nStates, nStates);
            t2p(i,j).zeros();

        }
    }
    cout << "Initialized t2 amplitudes." << endl;
}




