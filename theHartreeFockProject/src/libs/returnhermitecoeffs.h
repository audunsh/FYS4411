#ifndef RETURNHERMITECOEFFS_H
#define RETURNHERMITECOEFFS_H

#include <primitive.h>
#include <armadillo>

using namespace arma;

class ReturnHermiteCoeffs
{
public:
    ReturnHermiteCoeffs();
    field <cube> ReturnCoeffs(Primitive &Ga, Primitive &Gb);
    //field <mat> ReturnKineticMatrix();
    //double ReturnKineticIntegral();

private:
    field <cube> setup_E(const int &i_max, const  int &j_max, const int &k_max, const int &l_max, const int &m_max, const int &n_max);
    //void SetupKinteicIntegrals(const field<cube> &E, const double b);
    //double Sij(const field <cube> &E, const int ikj, const int ai, const int aj);
    void set_p(const double a,const double b);


    double p;               //  a+b
    //double Tab;             // Kinteic energy integral
    double pi = 4*atan(1);  // def pi
    //field <mat> T;          // Kinetic energy matrix


};

// Should use this for speed when the functionality si thouroughly validated!:
/*
inline const field<cube> &ReturnHermiteCoeffs::ReturnKineticIntegrals(){
    return T;
}
*/

#endif // ReturnHermiteCoeffs_H
