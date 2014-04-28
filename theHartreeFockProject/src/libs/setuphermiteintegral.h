#ifndef SETUPHERMITEINTEGRAL_H
#define SETUPHERMITEINTEGRAL_H

#include <primitive.h>
#include <boysfunction.h>
#include <armadillo>

using namespace arma;

class setupHermiteIntegral
{
public:
    setupHermiteIntegral(const Primitive &Ga, const Primitive &Gb, const vec corePos);
    field<cube> ReturnHermiteIntegral();

private:
    field <cube> Rtuv;
    void setupRtuv(const vec R, const vec F, const int &N, const int T, const int U, const int V);
    int n;
};

#endif // SETUPHERMITEINTEGRAL_H
