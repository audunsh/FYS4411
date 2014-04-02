#ifndef RETURNHERMITEPRODUCT_H
#define RETURNHERMITEPRODUCT_H

#include <primitive.h>
#include <armadillo>

using namespace arma;

class ReturnHermiteProduct
{
public:
    ReturnHermiteProduct(const Primitive &Ga, const Primitive &Gb );
private:
    void setup_E(field <cube> &E,
                 const int &i_max, const  int &j_max, const int &k_max, const int &l_max, const int &m_max, const int &n_max);
};

#endif // ReturnHermiteProduct_H
