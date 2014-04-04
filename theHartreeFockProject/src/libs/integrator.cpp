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


    vec P (3);                  // the middle point
    vec A (3);
    vec B (3);
    vec X_AB (3);
    vec X_PA (3);
    vec X_PB (3);
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
        X_AB[coord] = A[coord] - B[coord];
        X_PA[coord] = P[coord] - A[coord];
        X_PB[coord] = P[coord] - B[coord];
    }

    mu = a*b/p;

    /******************************************************************************/
    // setup Eij, Ekl, Emn for j = l = n = 0;


    double E_m,E_p;
    double Eij0,Ekl0,Emn0;
    /*******************************************************************************************
     *                          Eij. for the x-component:                                     */

    // for Eij
    if (i==0 && j==i) {
        Eij0 = exp(-mu*X_AB(0)*X_AB(0));
    }
    else {
        cube Eij (i+1,j+1,i+j);
        // initial values:
        //  i j t
        Eij(0,0,0) = exp(-mu*X_AB(0)*X_AB(0));

        if (i>0) {
            for (int ii = 0; ii < i; ++ii) {            // 0,1 // two values to find the third; i = 2 :-)
                for (int t = 0; t < i+j; ++t) {         // 0,1,2,3 // 4 values for i = j = 2.

                    if ((t-1) < 0) { E_m = 0;}
                    else { E_m = Eij(ii,0,t-1);}

                    if ((t+1) > ii) { E_p = 0;}
                    else { E_p = Eij(ii,0,t+1);}

                    if (t>ii) { Eij(ii,0,t) = 0.0;}

                Eij(ii+1,0,t) = 1/(2*p)*E_m + X_PA(0)*Eij(ii,0,t) + (t+1)*E_p;
                }
            }
        }
        // if i = 0 and j > 0 :
        else {
            for (int jj = 0; jj < j; ++jj) {
                for (int t = 0; t < i+j; ++t) {

                    if ((t-1) < 0) { E_m = 0;}
                    else {E_m = Eij(0,jj,t-1);}

                    if (t+1 > jj) { E_p = 0;}
                    else { E_p = Eij(0,jj,t+1);}

                    if (t > jj) { Eij(0,jj,t) = 0.0;}

                    Eij(0,jj+1,t) = 1/(2*p)*E_m + X_PB(0)*Eij(0,jj,t) + (t+1)*E_p;
                }
            }
        }

        // if i > 0 and j > 0.
        for (int ii = 1; ii <= i; ++ii) {
            for (int jj = 0; jj < j; ++jj) {
                for (int t = 0; t < i+j; ++t) {

                    if ((t-1) < 0) { E_m = 0;}
                    else {E_m = Eij(ii,jj,t-1);}

                    if (t+1 > ii+jj) { E_p = 0;}
                    else { E_p = Eij(ii,jj,t+1);}

                    if (t > ii+jj) { Eij(ii,jj,t) = 0.0;}

                    Eij(ii,jj+1,t) = 1/(2*p)*E_m + X_PB(0)*Eij(ii,jj,t) + (t+1)*E_p;
                }
            }
        }

        Eij0 = Eij(i,j,0);
    }

    /*******************************************************************************************
     *                             Ekl. for the y-component:                                  */

    if (k==0 && l==k) {
        Ekl0 = exp(-mu*X_AB(1)*X_AB(1));
    }
    else {
        Ekl0 = 1.0;
        cube Ekl (k+1,l+1,k+l);
        //  k l t
        Ekl(0,0,0) = exp(-mu*X_AB(1)*X_AB(1));


        if (k > 0) {
            for (int kk = 0; kk < k; ++kk) {
                for (int t = 0; t < k+l; ++t) {
                    if ((t-1) < 0) { E_m = 0;}
                    else { E_m = Ekl(kk,0,t-1);}

                    if (t+1 > kk) { E_p = 0;}
                    else { E_p = Ekl(kk,0,t+1);}

                    if (t>kk) { Ekl(kk,0,t) = 0.0;}

                    Ekl(kk+1,0,t) = 1/(2*p)*E_m + X_PA(1)*Ekl(kk,0,t) + (t+1)*E_p;
                }
            }
        }
        else {
            for (int ll = 0; ll < l; ++ll) {
                for (int t = 0; t < k+l; ++t) {

                    if ((t-1) < 0) { E_m = 0;}
                    else {E_m = Ekl(0,ll,t-1);}

                    if (t+1 > ll) { E_p = 0;}
                    else { E_p = Ekl(0,ll,t+1);}

                    if (t > ll) { Ekl(0,ll,t) = 0.0;}

                    Ekl(0,ll+1,t) = 1/(2*p)*E_m + X_PB(1)*Ekl(0,ll,t) + (t+1)*E_p;
                }
            }
        }
        for (int kk = 1; kk <= k; ++kk) {
            for (int ll = 0; ll < l; ++ll) {
                for (int t = 0; t < k+l; ++t) {

                    if ((t-1) < 0) { E_m = 0;}
                    else {E_m = Ekl(kk,ll,t-1);}

                    if (t+1 > kk+ll) { E_p = 0;}
                    else { E_p = Ekl(kk,ll,t+1);}

                    if (t > kk+ll) { Ekl(kk,ll,t) = 0.0;}

                    Ekl(kk,ll+1,t) = 1/(2*p)*E_m + X_PB(1)*Ekl(kk,ll,t) + (t+1)*E_p;
                }
            }
        }

        Ekl0 = Ekl(k,l,0);
    }

    /*******************************************************************************************
     *                         Emn. for the z-component:                                       */

    if (m==0 && n == m) {
        Emn0 = exp(-mu*X_AB(2)*X_AB(2));
    }
    else {
        cube Emn (m+1,n+1,m+n);
        //  m n t
        Emn(0,0,0) = exp(-mu*X_AB(2)*X_AB(2));

        if (m > 0) {
            for (int mm = 0; mm < m; ++mm) {
                for (int t = 0; t < m+n; ++t) {
                    if ((t-1) < 0) { E_m = 0;}
                    else {E_m = Emn(mm,0,t-1);}

                    if (t+1 > mm) {  E_p = 0;}
                    else { E_p = Emn(mm,0,t+1);}

                    if (t>mm) { Emn(mm,0,t) = 0.0;}

                    Emn(mm+1,0,t) = 1/(2*p)*E_m + X_PA(2)*Emn(mm,0,t) + (t+1)*E_p;
                }
            }
        }
        else {
            for (int nn = 0; nn < n; ++nn) {
                for (int t = 0; t < m+n; ++t) {

                    if ((t-1) < 0) { E_m = 0;}
                    else { E_m = Emn(0,nn,t-1);}

                    if (t+1 > nn) {E_p = 0;}
                    else { E_p = Emn(0,nn,t+1);}

                    if (t > nn) { Emn(0,nn,t) = 0.0;}

                    Emn(0,nn+1,t) = 1/(2*p)*E_m + X_PB(2)*Emn(0,nn,t) + (t+1)*E_p;
                }
            }
        }

        for (int mm = 1; mm <= m; ++mm) {
            for (int nn = 0; nn < n; ++nn) {
                for (int t = 0; t < m+n; ++t) {

                    if ((t-1) < 0) { E_m = 0;}
                    else { E_m = Emn(mm,nn,t-1);}

                    if (t+1 > mm+nn) {E_p = 0;}
                    else { E_p = Emn(mm,nn,t+1);}

                    if (t > mm+nn) { Emn(mm,nn,t) = 0.0;}

                    Emn(mm,nn+1,t) = 1/(2*p)*E_m + X_PB(2)*Emn(mm,nn,t) + (t+1)*E_p;
                }
            }
        }
        Emn0 = Emn(m,n,0);
    }

/*
    cout << "---------------------------------------------" << endl;
    cout << "i= " << i << " j= " << j << " k= " << k << " l= " << l << " m= " << m << " n= " << n << endl;


    for (int jj = 0; jj <= j; ++jj) {
        cout << "----------------------------------------------" << endl;
        cout << "-------------------  j=" << jj << " ---------------------" << endl;
        cout << "----------------------------------------------" << endl;
        for (int t = 0; t < i+j; ++t) {
            cout << "t=" << t  << " ";
            for (int ii = 0; ii <= i; ++ii) {
                cout << Eij(ii,jj,t) << " " ;
            }
            cout << endl;
        }
    }

    for (int ii = 0; ii<= i; ++ii) {
        cout << "----------------------------------------------" << endl;
        cout << "-------------------  i=" << ii << " ---------------------" << endl;
        cout << "----------------------------------------------" << endl;
        for (int t = 0; t < i+j; ++t) {
            cout << "t=" << t  << " ";
            for (int jj = 0; jj <= j; ++jj) {
                cout << Eij(ii,jj,t) << " " ;
            }
            cout << endl;
        }
    }


    // Have now found the cubic matrices Eij,Ekl,Emn.
    // it's time to find Sab = EijEklEmn(pi/p)^(3/2)
    cout << "---------------------------------------------------" << endl;
    cout << "Eij(i,j,0)= "<< Eij(i,j,0) << endl;
    cout << "Ekl(k,l,0)= "<< Ekl(k,l,0) << endl;
    cout << "Emn(m,n,0)= "<< Emn(m,n,0) << endl;
    cout << "X_AB=  " << X_AB(0) << " " << X_AB(1) << " " << X_AB(2) << endl;
    cout << "X_PA=  " << X_PA(0) << " " << X_PA(1) << " " << X_PA(2) << endl;
    cout << "X_PB=  " << X_PB(0) << " " << X_PB(1) << " " << X_PB(2) << endl;
    cout << "A=     " << A(0) << " " << A(1) << " " << A(2) << endl;
    cout << "B=     " << B(0) << " " << B(1) << " " << B(2) << endl;
*/

    cout << "Eij(i,j,0) = " << Eij0 << "  Eij(k,l,0) = " << Ekl0 << "  Emn(m,n,0) = " << Emn0 << endl;
    //double Sab = Eij(i,j,0)*Ekl(k,l,0)*Emn(m,n,0)*pow((pi/p),3.0/2);
    double Sab = Eij0*Ekl0*Emn0*pow((pi/p),3.0/2);
    return Sab;
}

void integrator::setupHermiteCoefficients(){
}