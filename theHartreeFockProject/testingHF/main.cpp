#include <unittest++/UnitTest++.h>
#include <myclass.h>
#include <lib.h>
#include <basis.h>
#include <boysfunction.h>
#include <integrator.h>
#include <primitive.h>
#include <hfsolve.h>

TEST(MyMath) {
  MyClass my;
  CHECK(my.addition(3,4) == 7);
}

TEST(HFSolve){
    // Test HFSolve using the known value of the groundstate energy of Beryllium
    // as a controll.

    //set number of protons and electrons
    int Z = 4;   // Number of protons
    int N = 4;   // Number of electrons
    int Ns = 4;  // 6 states

    basis Bs(3, 0);              //creating the basis object
    string filename;
    filename = "m_elements_c.dat";
    Bs.read("m_elements_c.dat", Z); //reading basis from file
    Bs.set_orthonormal(true);

    //Solving for N,Z with the provided basis
    HFSolve object (Z,N);
    double E = object.Solve(Bs);

    double E_GroundState_Beryllium = -14.6;
    double margin = 0.1;
    bool Bool = false;

    if (E_GroundState_Beryllium <= E){
        if (E - E_GroundState_Beryllium <= margin){
            Bool = true;
            cout << "HFSolve worked for groundstate of Beryllium" << endl;
            cout << "margin: " << E - E_GroundState_Beryllium << endl;
        }
    }
    if (Bool == false){
        cout << "diff energy of groundstate Beryllium: " << E - E_GroundState_Beryllium << endl;
        cout << "numerically E= " << E << " Ground State Beryllium E= " << E_GroundState_Beryllium << endl;
    }
    CHECK(Bool == true);
}


int main() {

    return UnitTest::RunAllTests();
}
