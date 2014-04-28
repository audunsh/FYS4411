#include "setuphermiteintegral.h"

// Rpc = P-C.
// p = a+b
// t = i+j, u = k+l, v = m+n
setupHermiteIntegral::setupHermiteIntegral(const double &p, const vec &Rpc, const int &t, const int &u, const int &v){
    Rtuv = field <mat> (3);
    n = t+u+v;
    double R = Rpc(0)*Rpc(0) + Rpc(1)*Rpc(1) +Rpc(2)*Rpc(2);
    double Fn, R0;
    BoysFunction F0(n);
    F0.setx(p*R*R);
    Fn = F0.returnValue(n);
    R0 = pow((-2.0*p),(double) n)*Fn;

    setupRtuv(Rpc,R0,n,t,u,v);
}

field <mat> setupHermiteIntegral::ReturnHermiteIntegral(){
    return Rtuv;
}

void setupHermiteIntegral::setupRtuv(const vec R, const double &R0, const int &N, const int T, const int U, const int V){
    // Doing the recurcion, and set up Rtuv.
    double X,Y,Z;

    Rtuv(0) = zeros <mat> (N+2,T+2);
    Rtuv(1) = zeros <mat> (N+2,U+2);
    Rtuv(2) = zeros <mat> (N+2,V+2);

    X = R(0);
    Y = R(1);
    Z = R(2);

    Rtuv(0)(N+1,1) = R0;

    cout << "R0=" << R0 << endl;
    cout << "n=" << (T+U+V) << " t=" << T << " u=" << U << " v=" << V << endl;
    cout << "--------------------R0-Rnt-----------------------" << endl;
    for (int n = N; n >= 0 ; --n) {
        for (int t = 1; t < T+1; ++t) {
            Rtuv(0)(n,t+1) = t*Rtuv(0)(n+1,t-1) + X*Rtuv(0)(n+1,t);
            cout << Rtuv(0)(n,t+1) << " " ;
        }
    }
    cout << endl;
    cout << "--------------------R0-Rnu-----------------------" << endl;
    for (int n = N; n >= 0; --n) {
        for (int u = 1; u < U+1; ++u) {
            Rtuv(1)(n,u+1) = u*Rtuv(1)(n+1,u-1) + Y*Rtuv(1)(n+1,u);
            cout << Rtuv(1)(n,u+1) << " " ;

        }
    }
    cout << endl;
    cout << "--------------------R0-Rnv-----------------------" << endl;
    for (int n = N; n >= 0; --n) {
        for (int v = 1; v < V+1; ++v) {
            Rtuv(2)(n,v+1) = v*Rtuv(2)(n+1,v-1) + Z*Rtuv(2)(n+1,v);
            cout << Rtuv(2)(n,v+1) << " " ;

        }
    }
    cout << endl;

}
