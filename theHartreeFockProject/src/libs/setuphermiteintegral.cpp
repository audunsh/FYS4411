#include "setuphermiteintegral.h"

// Rpc = P-C.
// p = a+b
// t = i+j, u = k+l, v = m+n
setupHermiteIntegral::setupHermiteIntegral(const double &p, const vec &Rpc, const vec F, const int &t, const int &u, const int &v){

    n = t+u+v;
    Rtuv = field <cube> (n+1);
    double R = Rpc(0)*Rpc(0) + Rpc(1)*Rpc(1) +Rpc(2)*Rpc(2);
    double Fn, R0;
    BoysFunction F0(n);
    F0.setx(p*R*R);
    Fn = F0.returnValue(n);
    R0 = pow((-2.0*p),(double) n)*Fn;

    setupRtuv(Rpc,R0,F,n,t,u,v);
}

field <cube> setupHermiteIntegral::ReturnHermiteIntegral(){
    return Rtuv;
}

void setupHermiteIntegral::setupRtuv(const vec R, const double &R0, const vec F, const int &N, const int T, const int U, const int V){
    // Doing the recurcion, and set up Rtuv.
    double X,Y,Z;

    for (int n = 0; n < N; ++n) {
        Rtuv(n) = zeros(T+1,U+1,V+1);
    }

    X = R(0);
    Y = R(1);
    Z = R(2);

    for (int n = 0; n < N; ++n) {
        Rtuv(n)(1,1,1) = F(n);
    }

    cout << "R0=" << R0 << endl;
    cout << "n=" << (T+U+V) << " t=" << T << " u=" << U << " v=" << V << endl;
    cout << "--------------------R0-Rnt-----------------------" << endl;
    for (int n = N; n >= 0 ; --n) {
        for (int t = 1; t < T+1; ++t) {
            Rtuv(n)(t+1,0,0) = t*Rtuv(n+1)(t-1,0,0) + X*Rtuv(n+1)(t,0,0);
            cout << Rtuv(n)(t+1,0,0) << " " ;
        }
    }

    cout << endl;
    cout << "--------------------R0-Rnu-----------------------" << endl;
    for (int n = N; n >= 0; --n) {
        for (int t = 0; t < T+1; ++t) {
            for (int u = 1; u < U+1; ++u) {
                Rtuv(n)(t,u+1,0) = u*Rtuv(n+1)(t,u-1,0) + Y*Rtuv(n+1)(t,u,0);
                cout << Rtuv(n)(t,u+1,0) << " " ;

            }
        }
    }
    cout << endl;
    cout << "--------------------R0-Rnv-----------------------" << endl;
    for (int n = N; n >= 0; --n) {
        for (int t = 0; t < T+1; ++t) {
            for (int u = 0; u < U+1; ++u) {
                for (int v = 1; v < V+1; ++v) {
                    Rtuv(n)(t,u,v+1) = v*Rtuv(n+1)(t,u,v-1) + Z*Rtuv(n+1)(t,u,v);
                    cout << Rtuv(n)(t,u,v+1) << " " ;
                }
            }
        }
    }
    cout << endl;

}

