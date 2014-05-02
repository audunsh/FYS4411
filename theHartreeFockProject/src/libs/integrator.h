#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <primitive.h>
#include <boysfunction.h>

using namespace std;
using namespace arma;

class integrator
{
public:
    integrator(Primitive &pA, Primitive &pB);
    double overlapIntegral(Primitive &pA, Primitive &pB);
    void setupHermiteCoefficients();
    void setupEij();
    void setupEcd();
    void setupRtuv(vec3 &nucleiPos);
    void setupRtau();

    //the integrals
    double overlap();
    double kinetic();
    double pNuclei();
    double pp(Primitive &pC, Primitive &pD);

private:
    //double pi = 4*atan(1);
    field <cube> Eij;
    field <cube> Ecd;
    field <cube> Rtuv;
    field <cube> Rtau;
    //vec P, pAijk, pBijk, A,B,Xab,Xpa,Xpb;
    vec3 P, pAijk, pBijk, pCijk, pDijk,A,B,C,D,Xab,Xcd,Xpa,Xpb,Xqc,Xqd,Rpc,Sijk,Tijk,Q, Rpq;
    double a,b,c,d,p,mu, Xab2,Xcd2,wA,wB,wC,wD,R, Rpc2,Rpq2,S,q, alpha;




};

#endif // INTEGRATOR_H
