#include "setuphermiteintegral.h"

// Rpc = P-C.
// p = a+b
// t = i+j, u = k+l, v = m+n
setupHermiteIntegral::setupHermiteIntegral(const double &p, const vec &Rpc, const int &t, const int &u, const int &v){
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
    Rtuv(0) = zeros(N+2,T+2);
    Rtuv(1) = zeros(N+2,U+2);
    Rtuv(2) = zeros(N+2,V+2);

    X,Y,Z = R(0), R(1), R(2);
    Rtuv(0)(N+1,1) = R0;

    for (int ni = 0; ni < N; ++ni) {
        n = N-ni;
        for (int t = 0; t < T+1; ++t) {
            Rtuv(0)(n,t+1) = t*Rtuv(0)(n+1,t-1) + X*Rtuv(0)(n+1,t);
        }
    }

    for (int ni = 0; ni < N; ++ni) {
         n = N-ni;
        for (int u = 0; u < U+1; ++u) {
            Rtuv(1)(n,u+1) = u*Rtuv(1)(n+1,u-1) + Y*Rtuv(1)(n+1,u);
        }
    }

    for (int ni = 0; ni < N; ++ni) {
        n = N-ni;
        for (int v = 0; v < V+1; ++v) {
            Rtuv(2)(n,v+1) = v*Rtuv(2)(n+1,v-1) + Z*Rtuv(2)(n+1,v);
        }
    }
}
