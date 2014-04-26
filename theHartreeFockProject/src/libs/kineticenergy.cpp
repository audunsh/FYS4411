#include "kineticenergy.h"

KineticEnergy::KineticEnergy(const field<cube> &E, const Primitive *Ga, const Primitive *Gb){

    m_E = E;
    A = Ga;
    B = Gb;
    b = Gb->exponent();
    p = Ga->exponent() + b;
}

double KineticEnergy::ReturnKineticIntegral(){
    int i,j,k,l,m,n;
    double Sij,Skl,Smn,Tij,Tkl,Tmn;
    double S = sqrt(pi/p);
    i = A->xExponent();
    k = A->yExponent();
    m = A->zExponent();

    j = B->xExponent();
    l = B->yExponent();
    n = B->zExponent();

    //cout << "i=" << i << " j=" << j <<" k=" << k <<" l=" << l << " m=" << m << " n=" << n << endl; // check; ok.

    Sij = m_E(0)(i,j,0);
    Skl = m_E(1)(k,l,0);
    Smn = m_E(2)(m,n,0);

    Tij = T_kin(S,i,j,0);
    Tkl = T_kin(S,k,l,1);
    Tmn = T_kin(S,m,n,2);

    cout << "----------------------------------------------------------------" << endl;
    cout << "i=" << i << " j=" << j <<" k=" << k <<" l=" << l << " m=" << m << " n=" << n << endl;
    cout << "Tij=" << Tij << " Tkl=" << Tkl << " Tmn=" << Tmn << endl;
    cout << "Sij=" << Sij << " Skl=" << Skl << " Smn=" << Smn << endl;

    return -0.5*(Tij*Skl*Smn + Sij*Tkl*Smn + Sij*Skl*Tmn);
}

double KineticEnergy::T_kin(const double S, const int i, const int j, const int cor){
    double Si_ = 0;
    double Si_p = 0;
    if ((j-2) >= 0) { Si_  = m_E(cor)(i,j-2,0);}
    if ((j+2) <= j) { Si_p = m_E(cor)(i,j+2,0);}

    double Tij = 4*b*b*Si_p - 2*b*(2*j+1)*m_E(cor)(i,j,0) + j*(j-1)*Si_;
    return Tij*S;
}
