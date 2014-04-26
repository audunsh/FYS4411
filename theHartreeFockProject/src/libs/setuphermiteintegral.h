#ifndef SETUPHERMITEINTEGRAL_H
#define SETUPHERMITEINTEGRAL_H

#include <armadillo>

using namespace arma;

class setupHermiteIntegral
{
public:
    setupHermiteIntegral(const double &p, const vec &Rpc);
    field <mat> ReturnHermiteIntegral();

private:
    field <mat> Rtuv;
};

#endif // SETUPHERMITEINTEGRAL_H
