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
    void setupRtuv(const double R0, const int n, const int t, const int u, const int v);
    int n;
};

#endif // SETUPHERMITEINTEGRAL_H
