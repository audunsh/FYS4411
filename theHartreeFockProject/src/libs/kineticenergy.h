#ifndef KINETICENERGY_H
#define KINETICENERGY_H

#include <armadillo>
#include <primitive.h>

using namespace arma;

class KineticEnergy
{
public:
    KineticEnergy(const field <cube> &E, const Primitive *Ga, const Primitive *Gb);
    double ReturnKineticIntegral();
private:
    double T_kin(const double S, const int i, const int j, const int cor);
    field <cube> m_E;
    const Primitive *A;
    const Primitive *B;
    const double pi = 4*atan(1);
    double p;
    double b;

};

#endif // KINETICENERGY_H
