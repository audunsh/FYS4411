#ifndef RETURNHERMITECOEFFS_H
#define RETURNHERMITECOEFFS_H

#include <primitive.h>
#include <armadillo>

using namespace arma;

class ReturnHermiteCoeffs
{
public:
    ReturnHermiteCoeffs();
    field<cube> ReturnCoeffs(Primitive &Ga, Primitive &Gb);
private:
    void setup_E(field <cube> &E,
                 const int &i_max, const  int &j_max, const int &k_max, const int &l_max, const int &m_max, const int &n_max);
};

#endif // ReturnHermiteCoeffs_H
