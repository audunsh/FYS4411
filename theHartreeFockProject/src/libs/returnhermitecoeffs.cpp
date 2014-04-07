#include <returnhermitecoeffs.h>


/*
 * Integrator 2 is a start on optimizing Integrator, or the setup of the Hermite coefficients to be more correctly! */

ReturnHermiteCoeffs::ReturnHermiteCoeffs(){

    // construct:
    T = field <mat> (3);       //  Holds the kinitic energies
}

field <cube> ReturnHermiteCoeffs::ReturnCoeffs(Primitive &Ga, Primitive &Gb){

    vec P (3);                  // the middle point
    vec A (3);
    vec B (3);
    vec K_AB (3);
    vec K_PA (3);
    vec K_PB (3);
    double a,b,mu;
    int i,j,k,l,m,n;

    a = Ga.exponent();          // exponential constant.
    A = Ga.nucleusPosition();   // nucleus Ga position
    i = Ga.xExponent();
    k = Ga.yExponent();
    m = Ga.zExponent();

    b = Gb.exponent();          // exponential constant.
    B = Gb.nucleusPosition();   // nucleus Gb position
    j = Gb.xExponent();
    l = Gb.yExponent();
    n = Gb.zExponent();

    set_p(a,b); // set global variable p.

    for (int coord = 0; coord < 3; ++coord) {
        P[coord] = (a*A[coord] + b*B[coord])/p;
        K_AB[coord] = A[coord] - B[coord];
        K_PA[coord] = P[coord] - A[coord];
        K_PB[coord] = P[coord] - B[coord];
    }

    mu = a*b/p;

    /***********************************************************************************************
     *                  Setup E part 2 :-) - We need to put in the actual values in E.
     *  calling setup_E(E i,j,k,l,m,n) just initialize the correct size of E, and fills it with
     *  values 666, to make it easier (hopefully ) to debug.                                       */

    field <cube> E = setup_E(i,j,k,l,m,n);
    for (int cor = 0; cor < E.n_elem; ++cor) {
        E(cor)(0,0,0) = exp(-mu*K_AB(cor)*K_AB(cor)); // initial value
    }

    // now we have set the initial value, it is time to iterate for the others,
    // beginning with one one direction, and then fill up the other.

    // E(0) = E(i,j,t); E(i+1,j,t) = (1/2p)*E(i-1,j,t) + K_PB(0)E(i,j,t) + (t+1)*E(i+1,j,t)
    int iBp,jBp;
    double E_m,E_p;
    for (int cor = 0; cor < E.n_elem; ++cor) { // loop for Eij, Ekl, Emn

        for (int iB = 1; iB < E(cor).n_rows; ++iB) {
            for (int t = 0; t < E(cor).n_slices; ++t) {

                iBp = iB - 1;   // previous value
                if ((t-1) < 0) { E_m = 0;}
                else { E_m = E(cor)(iBp,0,t-1);}

                if ((t+1) > iBp) { E_p = 0;}
                else { E_p = E(cor)(iBp,0,t+1);}

                if (t>iBp) { E(cor)(iBp,0,t) = 0.0;}

                E(cor)(iB,0,t) = 1/(2*p)*E_m + K_PA(cor)*E(cor)(iBp,0,t) + (t+1)*E_p;

            }
        }

        // filling up the other dimentions, that is, for j,l,n
        for (int jB = 1; jB < E(cor).n_cols; ++jB) {
            for (int iB = 0; iB < E(cor).n_rows; ++iB) {
                for (int t = 0; t < E(cor).n_slices; ++t) {

                    jBp = jB - 1; // previous value

                    if ((t-1) < 0) { E_m = 0;}
                    else {E_m = E(cor)(iB,jBp,t-1);}

                    if (t+1 > iB+jBp) { E_p = 0;}
                    else { E_p = E(cor)(iB,jBp,t+1);}

                    if (t > iB+jBp) { E(cor)(iB,jBp,t) = 0.0;}

                    E(cor)(iB,jB,t) = 1/(2*p)*E_m + K_PB(cor)*E(cor)(iB,jBp,t) + (t+1)*E_p;

                }
            }
        }
    }

    SetupKinteicIntegrals(E,b);
    return E;
}

field <mat> ReturnHermiteCoeffs::ReturnKineticMatrix(){
    return T;
}

double ReturnHermiteCoeffs::ReturnKineticIntegral(){
    return Tab;
}

void ReturnHermiteCoeffs::SetupKinteicIntegrals(const field<cube> &E, const double b){
    /*
     * Calculate the Kinetic integrals Tij,Tkl,Tmn and store them in a field of cubes to be collected from ReturnKinteicIntegrals */

    double Si_,Si_p;
    int i_max,j_max,k_max,l_max,m_max,n_max,iMAX,jMAX;

    i_max = E(0).n_rows - 1;
    j_max = E(0).n_cols - 1;
    k_max = E(1).n_rows - 1;
    l_max = E(1).n_cols - 1;
    m_max = E(2).n_rows - 1;
    n_max = E(2).n_cols - 1;

    for (int cor = 0; cor < 3; ++cor) {
        iMAX = E(cor).n_rows;
        jMAX = E(cor).n_cols;
        T(cor) = zeros <mat> (iMAX,jMAX);

        for (int i = 0; i < iMAX; ++i) {
            for (int j = 0; j < jMAX; ++j) {

                Si_ = 0;
                Si_p = 0;
                if ( (j-2) > 0) {
                    Si_ = Sij(E,cor,i,j-2);
                }
                if ((j+2) < j_max) {
                    Si_p = Sij(E,cor,i,j+2);
                }

                T(cor)(i,j) = 4*b*b*Si_p - 2*b*(2*i + 1)*Sij(E,cor,i,j) + j*(j-1)*Si_;
            }
        }
    }

    double Tij = T(0)(i_max, j_max);
    double Tkl = T(1)(k_max, l_max);
    double Tmn = T(2)(m_max, n_max);

    double Sij = E(0)(i_max,j_max,0);
    double Skl = E(1)(j_max,l_max,0);
    double Smn = E(2)(m_max,n_max,0);

    Tab = -0.5*(Tij*Skl*Smn + Sij*Tkl*Smn + Skl*Tmn);
}


double ReturnHermiteCoeffs::Sij(const field<cube> &E, const int xyz, const int i, const int j){
    /*
     * Calculate the overlap integral for index i,j between particle a and b.
     *         Sij = E(i,j,0)*exp((pi/p)^(3/2))                                */

    return E(xyz)(i,j,0)*pow(pi/p,1.0/2);
}



// Set up the right dimentionality of E, and fill it with numbers 666 - the number of the beast - because E is a fucking beast!!!!
field <cube> ReturnHermiteCoeffs::setup_E(const int &i_max, const int &j_max, const int &k_max, const int &l_max, const int &m_max, const int &n_max){

    // Set up a field of three cubes, each with the right dimentionality.

    field<cube> E(3);

    int t_max,u_max,v_max;
    t_max = i_max+j_max;
    u_max = k_max+l_max;
    v_max = m_max+n_max;

    //int E_size = E.size();
    E(0) = cube(i_max+1,j_max+1,t_max+1);
    E(1) = zeros <cube> (k_max+1,l_max+1,u_max+1);
    E(2) = zeros <cube> (m_max+1,n_max+1,v_max+1);

    for (int t = 0; t < t_max+1; ++t) {
        for (int i = 0; i < i_max+1; ++i) {
            for (int j = 0; j < j_max+1; ++j) {
                E(0)(i,j,t) = 666;
            }
        }
    }
    for (int u = 0; u < u_max+1; ++u) {
        for (int k = 0; k < k_max+1; ++k) {
            for (int l = 0; l < l_max+1; ++l) {
                E(1)(k,l,u) = 666;
            }
        }
    }
    for (int v = 0; v < v_max+1; ++v) {
        for (int m = 0; m < m_max+1; ++m) {
            for (int n = 0; n < n_max+1; ++n) {
                E(2)(m,n,v) = 666;
            }
        }
    }

    return E;
}



void ReturnHermiteCoeffs::set_p(const double a,const double b){
    p = a+b;
}



