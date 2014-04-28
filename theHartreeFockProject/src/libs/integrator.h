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
    void setupRtuv(vec &nucleiPos);
    void setupRtau(vec &nucleiPos, Primitive &pC, Primitive &pD);

    //the integrals
    double overlap();
    double kinetic();
    double pNuclei();
    double pp();

private:
    //double pi = 4*atan(1);
    field <cube> Eij;
    field <cube> Rtuv;
    field <cube> Rtau;
    //vec P, pAijk, pBijk, A,B,Xab,Xpa,Xpb;
    vec3 P, pAijk, pBijk, A,B,Xab,Xpa,Xpb;
    double a,b,p,mu, Xab2;




};

#endif // INTEGRATOR_H
