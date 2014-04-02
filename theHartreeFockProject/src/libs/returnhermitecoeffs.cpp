#include <returnhermitecoeffs.h>


/*
 * Integrator 2 is a start on optimizing Integrator, or the setup of the Hermite coefficients to be more correctly! */



ReturnHermiteCoeffs::ReturnHermiteCoeffs(field <cube> & E, const Primitive &Ga, const Primitive &Gb){

    vec P (3);                  // the middle point
    vec A (3);
    vec B (3);
    vec K_AB (3);
    vec K_PA (3);
    vec K_PB (3);
    double a,b,p,mu;
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

    p = a+b;
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

    setup_E(E,i,j,k,l,m,n);
    for (int cor = 0; cor < E.n_elem; ++cor) {
        E(cor)(0,0,0) = exp(-mu*K_AB(cor)*K_AB(cor)); // initial value
    }

    // now we have set the initial value, it is time to iterate for the others,
    // beginning with one one direction, and then fill up the other.

    // E(0) = E(i,j,t); E(i+1,j,t) = (1/2p)*E(i-1,j,t) + K_PB(0)E(i,j,t) + (t+1)*E(i+1,j,t)

    double E_m,E_p;
    for (int cor = 0; cor < E.n_elem; ++cor) { // loop for Eij, Ekl, Emn

        for (int iB = 1; iB < E(cor).n_rows; ++iB) {
            for (int t = 0; t < E(cor).n_slices; ++t) {

                if ((t-1) < 0) { E_m = 0;}
                else { E_m = E(cor)(iB,0,t-1);}

                if ((t+1) > iB) { E_p = 0;}
                else { E_p = E(cor)(iB,0,t+1);}

                if (t>iB) { E(cor)(iB,0,t) = 0.0;}

                E(cor)(iB+1,0,t) = 1/(2*p)*E_m + K_PA(cor)*E(cor)(iB,0,t) + (t+1)*E_p;

            }
        }

        // filling up the other dimentions, that is, for j,l,n
        for (int jB = 1; jB < E(cor).n_cols; ++jB) {
            for (int iB = 0; iB < E(cor).n_rows; ++iB) {
                for (int t = 0; t < E(cor).n_slices; ++t) {

                    if ((t-1) < 0) { E_m = 0;}
                    else {E_m = E(cor)(iB,jB,t-1);}

                    if (t+1 > iB+jB) { E_p = 0;}
                    else { E_p = E(cor)(iB,jB,t+1);}

                    if (t > iB+jB) { E(cor)(iB,jB,t) = 0.0;}

                    E(cor)(iB,jB+1,t) = 1/(2*p)*E_m + K_PB(cor)*E(cor)(iB,jB,t) + (t+1)*E_p;

                }
            }
        }
    }

}

void ReturnHermiteCoeffs::setup_E(field <cube> &E,
                          const int &i_max, const int &j_max, const int &k_max, const int &l_max, const int &m_max, const int &n_max){

    // Set up a field of three cubes, each with the right dimentionality.

    int t_max,u_max,v_max;
    t_max = i_max+j_max;
    u_max = k_max+l_max;
    v_max = m_max+n_max;

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

}



