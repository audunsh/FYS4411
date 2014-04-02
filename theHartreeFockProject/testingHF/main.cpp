#include <unittest++/UnitTest++.h>
#include <myclass.h>
#include <lib.h>
#include <basis.h>
#include <boysfunction.h>
#include <integrator.h>
#include <primitive.h>
#include <hfsolve.h>
#include <iomanip>


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

    cout << "----------------------------------------------------------------" << endl;
    cout << "------------------------ TEST HFSolve ----------------------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
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
    cout << "---------------------------END TEST HFSolve---------------------------------" << endl;
    CHECK(Bool == true);
}


TEST(integator_overlap_integrals_1){

    double a,b,weight;
    int i,j,k,l,m,n;
    vec A,B;

    // PrimitiveA
    a = 0.2;
    weight = 1;
    i = k = m = 0;
    A = {1.2, 2.3, 3.4};

    // PrimitiveB
    b = 0.3;
    weight = 1;
    j = l = n = 0;
    B = {-1.3,1.4,-2.4};

    Primitive primitiveA(weight,i,k,m,a,A);
    Primitive primitiveB(weight,j,l,n,b,B);

    integrator integral;
    double Sab = integral.overlapIntegral(primitiveA,primitiveB);
    cout << "----------------------------------------------------------------" << endl;
    cout << "----------------------- TEST integrator 1 ------------------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << "Sab = " << Sab << " and should be 0.1191723635809" << endl;
    cout << "---------------------- END TEST Integrator 1 ---------------------" << endl;
    CHECK_CLOSE(1.191723635809e-01,integral.overlapIntegral(primitiveA,primitiveB),1e-5);
}

TEST(integrator_overlap_integral_2){


    double a,b,weight;
    int i,j,k,l,m,n;
    vec A,B;

    // PrimitiveA:
    a = 0.2;
    weight = 1;
    i = m = 0;
    k = 1;
    A = {1.2, 2.3, 3.4};

    // PrimitiveB:
    b = 0.3;
    weight = 1;
    j = 0;
    l = n = 1;
    B = {-1.3,1.4,-2.4};

    Primitive primitiveA(weight,i,k,m,a,A);
    Primitive primitiveB(weight,j,l,n,b,B);
    integrator integral;

    double Sab = integral.overlapIntegral(primitiveA,primitiveB);
    cout << "----------------------------------------------------------------" << endl;
    cout << "------------------------- TEST Integrator 2 ----------------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << " Sab = " << Sab << " and should be 0.2227321941537" << endl;
    cout << "---------------------- END TEST Integrator 2 ---------------------" << endl;
    CHECK_CLOSE(2.227321941537e-01, Sab, 1e-5);
}

TEST(itegrator_overlap_integral_3){


    double a,b,weight;
    int i,j,k,l,m,n;
    vec A,B;

    // PrimitiveA:
    a = 0.2;
    weight = 1;
    i = m = 0;
    k = 2;
    A = {1.2, 2.3, 3.4};

    // PrimitiveB:
    b = 0.3;
    weight = 1;
    j = l = 1;
    n = 0;
    B = {-1.3,1.4,-2.4};

    Primitive primitiveA(weight,i,k,m,a,A);
    Primitive primitiveB(weight,j,l,n,b,B);
    integrator integral;
    double Sab = integral.overlapIntegral(primitiveA,primitiveB);
    cout << "----------------------------------------------------------------" << endl;
    cout << "------------------------ TEST Integrator 3 -----------------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << " Sab = " << Sab << " and should be -0.07329386373895" << endl;
    cout << "---------------------- END TEST Integrator 3 -------------------" << endl;
    CHECK_CLOSE(-7.329386373895e-02,Sab, 1e-5);
}




int main() {

    return UnitTest::RunAllTests();
}
