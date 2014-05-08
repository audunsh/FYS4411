#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include <string>
#include <boysfunction.h>
#include <basis.h>
#include <integrator.h>
#include <hfsolve.h>
#include <contracted.h>

double pi = 4*atan(1);

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {
    /*
    //set number of protons and electrons
    int Z = 4;   // Number of protons
    int N = 4;   // Number of electrons
    int Ns = 4;  // 6 states

    basis BS(N, 0);              //creating the basis object
    string filename;
    filename = "m_elements_c.dat";
    BS.read(filename, Z); //reading basis from file
    BS.set_orthonormal(true);

    cout << "Energy of the ground state= " << E << endl;
    */
    basis BS(3); //set up a basis containing 3 contracted/orbitals
    BS.init_STO_3G("Be"); //initialize the STO-3G basis for the Beryllium atom
    BS.init_integrals();  //set up and solve the needed integrals to calculate overlapmatrix, single-particle interaction and two-particle interaction
    //output the result



    for(int p=0;p<3;p++){
        for(int q=0;q<3;q++){
            for(int r=0;r<3;r++){
                for(int s=0;s<3;s++){
                    cout << p << q << r << s << " " << BS.v(p,q)(r,s) << endl;
                }
            }
            cout << "    " << p << q << BS.h(p,q) << " " << BS.S(p,q);
        }
    }

    //HFSolve object (Z,N);
    //double E = object.Solve(BS);
    //cout << "Energy of the ground state= " << E << endl;
    return 0;
} // End: function output()
