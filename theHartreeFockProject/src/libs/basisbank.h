#ifndef BASISBANK_H
#define BASISBANK_H
#include <armadillo>
#include <basis.h>
#include <primitive.h>
using namespace std;
using namespace arma;

class basisbank{
public:
    basisbank(basis BS);
    basis bs;
    void add_STO6G_He(vec3 corePos);
    void add_STO6G_H(vec3 corePos);
    void add_STO6G_F(vec3 corePos);
    void add_6311G_H(vec3 corePos);
    void add_STO6G_B(vec3 corePos);
    void add_STO6G_Mg(vec3 corePos);
    void add_STO6G_Be(vec3 corePos);
    void add_STO6G_C(vec3 corePos);
    void add_STO6G_Na(vec3 corePos);
    void add_631ppG_H(vec3 corePos);
    void add_STO6G_O(vec3 corePos);
    void add_STO6G_N(vec3 corePos);
    void add_STO6G_Ne(vec3 corePos);
    void add_STO6G_Li(vec3 corePos);
};
#endif // BASISBANK_H
