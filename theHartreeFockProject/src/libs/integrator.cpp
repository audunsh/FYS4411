#include <integrator.h>
#include <primitive.h>
#include <armadillo>

using namespace std;
using namespace arma;

/***********************************************************************************
 *  Calculate integrals between primitive objects collected from class Primitive.
 *
 **********************************************************************************/
integrator::integrator(){

}

double integrator::overlapIntegral(Primitive Ga, Primitive Gb){

    // 1) setup Ga = Gikm(a,rA) and Gb = Gjln(b,rB).

    // 2) calculate the coefficients Eij, Ekl, Emn for the 0th recurcion

    // Sab = EijEklEmn(pi/p)^(3/2)  // this is the overlap integral between a and b. p = a+b


    // Eij = exp(-(ab/(a+b))X^2), X = X_{AB}

    //mA = a.weight();
    //mB = b.weight();
    //a.xExponent();
    //a.yExponent();

    //Eij = K_AB()

    vec P (3);
    vec A (3);
    vec B (3);
    vec X_AB (3);
    vec X_PA (3);
    vec X_PB (3);
    double a,b,p;
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
        X_AB[coord] = A[coord] - B[coord];
        X_PA[coord] = P[coord] - A[coord];
        X_PB[coord] = P[coord] - B[coord];
    }



    /******************************************************************************/
    // setup Eij, Ekl, Emn for j = l = n = 0;

    //    Ers(e,f,t)
    cube Eij (i,j,i+j);
    cube Ekl (k,l,k+l);
    cube Emn (m,n,m+n);

    // initial values:
    //    t  i  j
    Eij(0,0,0) = X_AB(0);
    //    t  k  l
    Ekl(0,0,0) = X_AB(1);
    //    t  m  n
    Emn(0,0,0) = X_AB(2);

    cout << "---------------------------------------------" << endl;
    cout << "i= " << i << " j= " << j << " k= " << k << " l= " << l << " m= " << m << " n= " << n << endl;

    for (int kapp = 0; kapp < i+j; ++kapp) {
        for (int kupp = 0; kupp < i; ++kupp) {
            for (int pupp = 0; pupp < j; ++pupp) {
                cout << Eij(kupp,pupp,kapp) << " " ;
            }
            cout << endl;
        }
        cout << "----------------------------------------------------" << endl;
        cout << "              NEW t VALUE STARTS HERE               " << endl;
        cout << "----------------------------------------------------" << endl;
    }

    // for Eij
    double E_m,E_p;
    for (int t = 0; t < i+j; ++t) {
        for (int ii = 0; ii < i; ++ii) {
            if ((t-1) < 0) {
                E_m = 0;
            }
            if (t+1 > ii) {
                E_p = 0;
            }
            else {
                E_m = Eij(ii,0,t-1);
                E_p = Eij(ii,0,t+1);
            }
            Eij(ii+1,0,t) = 1/(2*p)*E_m + X_PA(0)*Eij(ii,0,t) + (t+1)*E_p;
        }
    }

    for (int t = 0; t < k+l; ++t) {
        for (int kk = 0; kk < k; ++kk) {
            if ((t-1) < 0) {
                E_m = 0;
            }
            if (t+1 > kk) {
                E_p = 0;
            }
            else {
                E_m = Ekl(t-1,kk,0);
                E_p = Ekl(t+1,kk,0);
            }
            Ekl(t,kk+1,0) = 1/(2*p)*E_m + X_PA(1)*Ekl(t,kk,0) + (t+1)*E_p;
        }
    }

    for (int t = 0; t < m+n; ++t) {
        for (int mm = 0; mm < m; ++mm) {
            if ((t-1) < 0) {
                E_m = 0;
            }
            if (t+1 > mm) {
                E_p = 0;
            }
            else {
                E_m = Emn(t-1,mm,0);
                E_p = Emn(t+1,mm,0);
            }
            Emn(t,mm+1,0) = 1/(2*p)*E_m + X_PA(2)*Emn(t,mm,0) + (t+1)*E_p;
        }
    }


    // can now iterate for j,l,n because we now know Ers(:,0,:)

    for (int ii = 0; ii < i; ++ii) {
        for (int t = 0; t < i+j; ++t) {
            for (int jj = 1; jj < j; ++jj) {  // We already know the values for i=j=0 for all t, so we start at j=1.

                if ((t-1) < 0) {
                    E_m = 0;
                }
                if (t+1 > ii+jj) {
                    E_p = 0;
                }
                else {
                    E_m = Eij(t-1,ii,jj);
                    E_p = Eij(t+1,ii,jj);
                }
                Eij(t,ii,jj+1) = 1/(2*p)*E_m + X_PB(0)*Ekl(t,ii,jj) + (t+1)*E_p;
            }
            }
        }

    for (int kk = 0; kk < k; ++kk) {
        for (int t = 0; t < i+j; ++t) {
            for (int ll = 1; ll < l; ++ll) {

                if ((t-1) < 0) {
                    E_m = 0;
                }
                if (t+1 > kk+ll) {
                    E_p = 0;
                }
                else {
                    E_m = Ekl(t-1,kk,ll);
                    E_p = Ekl(t+1,kk,ll);
                }
                Eij(t,kk,ll+1) = 1/(2*p)*E_m + X_PB(1)*Ekl(t,kk,ll) + (t+1)*E_p;
            }
        }
    }


    for (int mm = 0; mm < m; ++mm) {
        for (int t = 0; t < i+j; ++t) {
            for (int nn = 1; nn < n; ++nn) {

                if ((t-1) < 0) {
                    E_m = 0;
                }
                if (t+1 > mm+nn) {
                    E_p = 0;
                }
                else {
                    E_m = Eij(t-1,mm,nn);
                    E_p = Eij(t+1,mm,nn);
                }
                Eij(t,mm,nn+1) = 1/(2*p)*E_m + X_PB(2)*Ekl(t,mm,nn) + (t+1)*E_p;
            }
        }
    }


    // Have now found the cubic matrices Eij,Ekl,Emn.
    // it's time to find Sab = EijEklEmn(pi/p)^(3/2)

    double Sab = Eij(0,i,j)*Ekl(0,k,l)*Emn(0,m,n)*pow((pi/p),3.0/2);
    return Sab;
}

void integrator::setupHermiteCoefficients(){
}
