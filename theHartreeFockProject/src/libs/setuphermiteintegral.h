#ifndef SETUPHERMITEINTEGRAL_H
#define SETUPHERMITEINTEGRAL_H

#include <boysfunction.h>
#include <armadillo>

using namespace arma;

class setupHermiteIntegral
{
public:
    setupHermiteIntegral(const double &p, const vec &Rpc, const vec F, const int &t, const int &u, const int &v);
    field<cube> ReturnHermiteIntegral();

private:
    field <cube> Rtuv;
    void setupRtuv(const vec R, const double &R0, const vec F, const int &N, const int T, const int U, const int V);
    int n;
};

#endif // SETUPHERMITEINTEGRAL_H
