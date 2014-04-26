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
<<<<<<< HEAD
    //double Tab;             // Kinteic energy integral
    double pi = 4*atan(1);  // def pi
    //field <mat> T;          // Kinetic energy matrix
=======
    double Tab;             // Kinteic energy integral
    double pi;  // def pi
    field <mat> T;          // Kinetic energy matrix
>>>>>>> ba732e7b13827689a65e8ed7735e8c3b4eac1940

};

// Should use this for speed when the functionality si thouroughly validated!:
/*
inline const field<cube> &ReturnHermiteCoeffs::ReturnKineticIntegrals(){
    return T;
}
*/

#endif // ReturnHermiteCoeffs_H
