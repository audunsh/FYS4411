/*
Calculations on the Hydrogen atom
- using the radial part of the Shroedinger equation to calculate the eigenvalues of the matrix equation
*/

#include <iostream>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include <lib.h>

using namespace std;
using namespace arma;

#define pi 4*arctan(1);          // pi
#define hbar 6.625*pow(10,-34)   // Planck's reduced constant
#define masse 9.109*pow(10,-31)  // mass of electron

int comp(double, double);
void output(double, double, int, double *);

int main(int argc, char* argv[]) {

    double k, m, Rmax, Rmin, alpha, h, lambda, eu;
    int l, Z, Nstep;
    k = 1.44;        // eVnm
    m = masse;       // mass of electron
    Rmin = 0;        // min radius (Bohr distances)
    Rmax = 10.0;     // max radius (Bohr distances)
    l = 0;           //
    Z = 1;           // atomic number
    Nstep = 100;      // number of intervalls
    alpha = hbar*hbar/(k*Z*m);
    h = (Rmax-Rmin)/Nstep;
    eu = -1/(2*h*h);


    //rowvec x(Nstep),d(Nstep),sub_d(Nstep);   // vectors containing diagonal and off-diagonal elements
    //mat U = zeros(Nstep,Nstep); // matrix containing the normalized eigenvectors on the k'th row, correstponding to d[k]

    double *x, *d, *sub_d, **U;
    d = new double[Nstep-2];
    x = new double[Nstep-2];
    sub_d = new double[Nstep-2];
    U = new double*[Nstep-2];

    for (int i=0;i<(Nstep-2);i++){
        U[i] = new double[Nstep-2];
    }
    for (int i=0;i<(Nstep-2);i++){
        // does not include the boundaries!
        x[i] = Rmin + (i+1)*h;
        sub_d[i] = eu;
        d[i] = (1/(h*h) + (l*(l+1))/(2*x[i]*x[i]) - 1/x[i]);
        U[i][i] = 1.0;
    }


    // Use the library function tqli to solve the eigenvalue problem
    tqli(d,sub_d,Nstep-2,U);

    /*
    // have a look at the results:
    for (int i=0;i<(Nstep-2);i++){
        for (int j=0;j<(Nstep-2);j++){
            cout << U[i][j] << " ";
        }
        cout << endl;
    }
    */
    cout << "-----------------------------------------------------------" << endl;
    cout << "x"<< endl;
    for (int i = 0; i < (Nstep-2); ++i) {
        cout << x[i] <<  endl;
    }
    cout << "-----------------------------------------------------------" << endl;
    cout << "sub_d"<< endl;
    for (int i = 0; i < (Nstep-2); ++i) {
        cout << sub_d[i] << endl;
    }
    cout << "-----------------------------------------------------------" << endl;
    cout << "d"<< endl;
    for (int i = 0; i < (Nstep-2); ++i) {
        cout << d[i] << endl;
    }
    cout << "-----------------------------------------------------------" << endl;

    //find the lowest eigenvalue in the d-array:
    double min_eig = d[0];
    for (int i = 1; i < Nstep - 2 ; ++i) {
        if (d[i] < d[i-1]){
            if (d[i] < d[i+1]){
                if (d[i] < min_eig){
                    min_eig = d[i];
                }
            }
        }
    }
    cout << min_eig << endl;

    // this (lower) part is not finished yet, I do not know why it isn't working, but the debugger say
    // that the reference to the functions comp() and to output() is undefined.

    /*Sort the matrix values, need to store the eigenvalues from smallest to largest in a vector
      where the indexes of the original positions of the eigenvalues is the indexes of the columns
      in the unsorted matrix containing the corresponding eigenvectors.
      */

    qsort(d,(UL) Nstep - 1,sizeof(double),(int(*)(const void *,const void *))comp);

    // write to outputfile
    output(Rmin,Rmax,Nstep,*d);

    // free memory:
    free_matrix((void **) U);
    delete [] x;
    delete [] d;
    delete [] sub_d;

    return 0;
} // End: main()




int comp(const double *val_1, const double *val_2){
  if((*val_1) <= (*val_2))       return -1;
  else  if((*val_1) > (*val_2))  return +1;
  else                     return  0;
} // End: function comp()

void output(double r_min,double r_max, int N_step, double *d){
    ofstream myfile;
    myfile.open("results.txt");
    myfile << "RESULTS:" << endl;
    myfile << setiosflags(ios::showpoint | ios::uppercase);
    myfile << "Rmin = " << setw(15) << setprecision(8) << r_min << endl;
    myfile << "Rmax = " << setw(15) << setprecision(8) << r_max << endl;
    myfile << "Number of steps = " << setw(15) << N_step << endl;
    myfile << "Eigenvalues: " << endl;
    for (int i=0;i<N_step-1;i++){
        myfile << setw(15) << setprecision(8) << d[i] << endl;
    }
    myfile.close();
} // End: function output()
    return 0;
}









