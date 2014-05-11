#include <unittest++/UnitTest++.h>
#include <lib.h>
#include <basis.h>
#include <boysfunction.h>
#include <integrator.h>
#include <primitive.h>
#include <hfsolve.h>
#include <iomanip>


double pi = 4*atan(1);

TEST(HFSolve){

    // Test HFSolve using the known value of the groundstate energy of Beryllium
    // as a controll.

    //set number of protons and electrons
    int Z = 4;   // Number of protons
    int N = 4;   // Number of electrons
    int Ns = 4;  // 6 states

    basis Bs;              //creating the basis object
    string filename;
    filename = "/home/goranbs/goran/CompPhys/FYS4411\ -\ CompPhys2/build-theHartreeFockProject-Desktop-Release/testingHF/m_elements_c.dat";
    Bs.read(filename, Z); //reading basis from file
    Bs.set_orthonormal();

    //Solving for N,Z with the provided basis
    HFSolve object (Z,N);
    double E = object.Solve(Bs);

    double E_GroundState_Beryllium = -14.6;
    double margin = 0.1;
    bool Bool = false;

    cout << "----------------------------------------------------------------" << endl;
    cout << "------------------------ TEST HFSolve --------------------------" << endl;
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
    cout << "________________________________________________________________" << endl;

    CHECK(Bool == true);
}


TEST(integrator1){
    BoysFunction boys(2);

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

    integrator AB(primitiveA,primitiveB,boys);

    double Sab = AB.overlap(); //Eab(0)(i,j,0)*Eab(1)(k,l,0)*Eab(2)(m,n,0)*pow(pi/(a+b),3.0/2);
    cout << "----------------------------------------------------------------" << endl;
    cout << "----------------------- TEST integrator 1 ----------------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << " Sab = " << Sab << " and should be:" << endl;
    cout << " Sab = 0.1191723635809" << endl;
    cout << "________________________________________________________________" << endl;

    CHECK_CLOSE(1.191723635809e-01,Sab,1e-5);
}

TEST(integrator2){

    BoysFunction boys(2);

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
    integrator AB(primitiveA,primitiveB,boys);

    double Sab = AB.overlap();
    cout << "----------------------------------------------------------------" << endl;
    cout << "------------------------ TEST Integrator 2 ---------------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << " Sab = " << Sab << " and should be:" << endl;
    cout << " Sab = 0.2227321941537" << endl;
    cout << "________________________________________________________________" << endl;

    CHECK_CLOSE(2.227321941537e-01, Sab, 1e-5);
}

TEST(integrator3){
    BoysFunction boys(2);


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
    integrator AB(primitiveA,primitiveB,boys);

    double Sab = AB.overlap();
    cout << "----------------------------------------------------------------" << endl;
    cout << "----------------------- TEST Integrator 3 ----------------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << " Sab = " << Sab << " and should be:" << endl;
    cout << " Sab = -0.07329386373895" << endl;
    cout << "________________________________________________________________" << endl;

    CHECK_CLOSE(-7.329386373895e-02,Sab, 1e-5);
}

TEST(Kinetic_integral_1){
    BoysFunction boys(2);
    double a,b,weight;
    int i,j,k,l,m,n;
    vec A,B;

    // PrimitiveA
    a = 0.2;
    weight = 1;
    i = k = m = 0; // 0
    A = { 1.2, 2.3, 3.4 } ;

    // PrimitiveB
    b = 0.3;
    weight = 1;
    j = l = n = 0; // 0
    B = {-1.3, 1.4, -2.4 };

    Primitive primitiveA(weight,i,k,m,a,A);
    Primitive primitiveB(weight,j,l,n,b,B);

    integrator AB(primitiveA,primitiveB,boys);

    //double Tab = Coeffs.ReturnKineticIntegral();

    double Tab = AB.kinetic();
    cout << "----------------------------------------------------------------" << endl;
    cout << "                   TESTING THE KINTETIC INTEGRALS               " << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << "--------------------- TEST Kinetic Integral 1 ------------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << " Tab = " << Tab << " and should be:" << endl;
    cout << " Tab = -0.09678702680582" << endl;
    cout << "________________________________________________________________" << endl;

    CHECK_CLOSE( -9.678702680582e-02, Tab, 1e-5);
}

TEST(Kinetic_integral_2){
    BoysFunction boys(2);
    double a,b,weight;
    int i,j,k,l,m,n;
    vec A,B;

    // PrimitiveA:
    a = 0.2;
    weight = 1;
    i = m = 0; // 0
    k = 1;  // 1
    A = { 1.2, 2.3, 3.4 };

    // PrimitiveB:
    b = 0.3;
    weight = 1;
    j = 0; // 0
    l = n = 1; // 1
    B = {-1.3, 1.4, -2.4 };
    Primitive primitiveA(weight,i,k,m,a,A);
    Primitive primitiveB(weight,j,l,n,b,B);
    integrator AB(primitiveA,primitiveB,boys);
    double Tab = AB.kinetic();

    cout << "----------------------------------------------------------------" << endl;
    cout << "------------------- TEST Kinetic Integral 2 --------------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << " Tab = " << Tab << " and should be:" << endl;
    cout << " Tab = -0.08688217105502" << endl;
    cout << "________________________________________________________________" << endl;


    CHECK_CLOSE( -8.688217105502e-02, Tab, 1e-5);

}

TEST(Kinetic_integral_3){
    BoysFunction boys(2);
    double a,b,weight;
    int i,j,k,l,m,n;
    vec A,B;

    // PrimitiveA:
    a = 0.2;
    weight = 1;
    i = m = 0; // 0
    k = 2; // 2
    A = { 1.2, 2.3, 3.4 };

    // PrimitiveB:
    b = 0.3;
    weight = 1;
    j = l = 1; // 1
    n = 0; // 0
    B = {-1.3, 1.4, -2.4 };
    Primitive primitiveA(weight,i,k,m,a,A);
    Primitive primitiveB(weight,j,l,n,b,B);
    integrator AB(primitiveA,primitiveB,boys);
    double Tab = AB.kinetic();
    cout << "----------------------------------------------------------------" << endl;
    cout << "-------------------- TEST Kinetic Integral 3 -------------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << " Tab = " << Tab << " and should be:" << endl;
    cout << " Tab = -0.01598401092187" << endl;
    cout << "________________________________________________________________" << endl;

    CHECK_CLOSE( -1.598401092187e-02, Tab, 1e-5);
}




int main() {

    return UnitTest::RunAllTests();
}
