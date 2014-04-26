#ifndef SETUPHERMITEINTEGRAL_H
#define SETUPHERMITEINTEGRAL_H

#include <boysfunction.h>
#include <armadillo>

using namespace arma;

class setupHermiteIntegral
{
public:
    setupHermiteIntegral(const double &p, const vec &Rpc, const int &t, const int &u, const int &v);
    field <mat> ReturnHermiteIntegral();

private:
    field <mat> Rtuv;
    void setupRtuv(const vec R, const double &R0, const int &N, const int T, const int U, const int V);
    int n;
};

#endif // SETUPHERMITEINTEGRAL_H
