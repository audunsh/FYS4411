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

    setupRtuv(R0,n,t,u,v);
}

field <mat> setupHermiteIntegral::ReturnHermiteIntegral(){

    return Rtuv;
}

void setupHermiteIntegral::setupRtuv(const double R0, const int n, const int t, const int u, const int v){
    // Doing the recurcion, and set up Rtuv.
}
