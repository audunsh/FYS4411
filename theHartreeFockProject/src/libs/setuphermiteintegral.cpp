#include "setuphermiteintegral.h"

// Rpc = P-C.
// p = a+b
// t = i+j, u = k+l, v = m+n
setupHermiteIntegral::setupHermiteIntegral(const Primitive &Ga, const Primitive &Gb, const vec corePos){

    double p,R,R0;
    vec P,Rpc;
    int n,t,u,v;

    p = Ga.exponent() + Gb.exponent();
    P = (Ga.exponent()*Ga.nucleusPosition() + Gb.exponent()*Gb.nucleusPosition())/p;
    Rpc = P - corePos;

    t = Ga.xExponent() + Gb.xExponent();
    u = Ga.yExponent() + Gb.yExponent();
    v = Ga.zExponent() + Gb.zExponent();
    n = t+u+v+2;

    Rtuv = field <cube> (n);
    R = Rpc(0)*Rpc(0) + Rpc(1)*Rpc(1) +Rpc(2)*Rpc(2);

    BoysFunction boys(n-1);    // Her maa vi vere obs! Det funcker for n < 6 !!! Vet ikke hvorfor....
    boys.setx(p*R*R);

    vec F(n);
    for (int i = 0; i < n-2; ++i) {   // skal denne begynne paa null?
        F(i) = pow((-2.0*p),(double) i)*boys.returnValue(i);
    }

    setupRtuv(Rpc,F,n,t,u,v);
}

field <cube> setupHermiteIntegral::ReturnHermiteIntegral(){
    return Rtuv;
}

void setupHermiteIntegral::setupRtuv(const vec R, const vec F, const int &N, const int T, const int U, const int V){
    // Doing the recurcion, and set up Rtuv.
    double X,Y,Z;

    for (int n = 0; n < N; ++n) {
        Rtuv(n) = zeros(T+2,U+2,V+2);
    }

    X = R(0);
    Y = R(1);
    Z = R(2);

    for (int n = 0; n < N; ++n) {
        //cout << F(n) << " " ;
        Rtuv(n)(1,1,1) = F(n);
    }
    cout << "n=" << (T+U+V) << " t=" << T << " u=" << U << " v=" << V << endl;
    cout << "--------------------R0-Rnt-------------------------------------------" << endl;
    for (int n = N-2; n > 0 ; --n) {
        cout << "n= " << n << " : " ;
        for (int t = 1; t < T+1; ++t) {
            Rtuv(n)(t+1,1,1) = t*Rtuv(n+1)(t-1,1,1) + X*Rtuv(n+1)(t,1,1);
            cout << Rtuv(n)(t-1,1,1) << " " ;
        }
        cout << Rtuv(n)(T,1,1) << " " << endl;
    }


    cout << "--------------------R0-Rnu-------------------------------------------" << endl;
    for (int n = N-2; n > 0; --n) {
        cout << "n= " << n << " : " << endl;
        for (int t = 0; t < T+1; ++t) {
            cout << " t= " << t << " : " ;
            for (int u = 1; u < U+1; ++u) {
                Rtuv(n)(t,u+1,1) = u*Rtuv(n+1)(t,u-1,1) + Y*Rtuv(n+1)(t,u,1);
                cout << Rtuv(n)(t,u-1,1) << " " ;

            }
            cout << Rtuv(n)(t,U,1) << endl;
        }
        cout << endl;
    }
    cout << endl;
    cout << "--------------------R0-Rnv-------------------------------------------" << endl;
    for (int n = N-2; n > 0; --n) {
        cout << "n= " << n << " : " << endl;
        for (int t = 0; t < T+1; ++t) {
            cout << " t= " << t << " : " << endl;
            for (int u = 0; u < U+1; ++u) {
                cout << "   u= " << u << " : ";
                for (int v = 1; v < V+1; ++v) {
                    Rtuv(n)(t,u,v+1) = v*Rtuv(n+1)(t,u,v-1) + Z*Rtuv(n+1)(t,u,v);
                    cout << Rtuv(n)(t,u,v-1) << " " ;
                }
                cout << Rtuv(n)(t,u,V) << endl;
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;

}

