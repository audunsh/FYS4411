#include <unittest++/UnitTest++.h>
#include <lib.h>
#include <basis.h>
#include <boysfunction.h>
#include <integrator.h>
#include <primitive.h>
#include <hfsolve.h>
#include <iomanip>
#include <hartreefocksolver.h>


double pi = 4*atan(1);

TEST(H2_STO3G){
    //Test H2 hydrogen molecule to given precision using STO-3G basis
    basis BS;
    int nElectrons = 2;
    double nProtons = 2;
    BS.init_H2({2,1.3,0},{2,2.7,0}); //insert parameter dist here (calculation is however still off for molecules)
    BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction
    hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
    double E = object.solve();                          //solve for the given basis
    cout << "H2:" << E << endl;
    //cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (Approx. " << 27.212*E << " eV)" << endl;
    bool B = false;
    double margin = 0.1;
    double target = -1.175;
    if(abs(E-target)<=margin){
        B = true;
    }
    CHECK(B);
}

TEST(He_STO3G){
    //Test Helium atom to given precision using STO-3G basis
    basis BS;
    int nElectrons = 2;
    double nProtons = 2;
    //BS.init_STO_3G("He", nProtons);
    BS.init_He({0,0,0});
    BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction
    hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
    double E = object.solve();                          //solve for the given basis
    //cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (Approx. " << 27.212*E << " eV)" << endl;

    bool B = false;
    double margin = 0.1;
    double target = -2.904;
    if(abs(E-target)<=margin){
        B = true;
    }
    CHECK(B);
}

/*
TEST(O_STO3G){
    //Test Helium atom to given precision using STO-3G basis
    basis BS;
    int nElectrons = 8;
    double nProtons = 8;
    BS.init_STO_3G("O", nProtons);
    BS.init_integrals();  //set up and solve the needed integrals to calculate overlap matrix, single-particle interaction and two-particle interaction
    hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
    double E = object.solve();                          //solve for the given basis
    //cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (Approx. " << 27.212*E << " eV)" << endl;
    bool B = false;
    double margin = 0.1;
    double target = -2.904;
    if(abs(E-target)<=margin){
        B = true;
    }
    CHECK(B);
}

*/

TEST(Be_Hydrogenlike){
    //Test Beryllium atom using hydrogenlike basis (provided as a table)
    basis BS;
    int nElectrons = 4;
    double nProtons = 4;
    BS.init_HTO4(nProtons);
    hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
    double E = object.solve();                          //solve for the given basis
    //cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (Approx. " << 27.212*E << " eV)" << endl;
    bool B = false;
    double margin = 0.2;
    double target = -14.67;
    if(abs(E-target)<=margin){
        B = true;
    }
    CHECK(B);
}

TEST(He_Hydrogenlike){
    //Test Helium atom using hydrogenlike basis (provided as a table)
    basis BS;
    int nElectrons = 2;
    double nProtons = 2;
    BS.init_HTO4(nProtons);
    hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
    double E = object.solve();                          //solve for the given basis
    //cout << setprecision(10) << "Ground state energy:" << E << " atomic units. (Approx. " << 27.212*E << " eV)" << endl;
    bool B = false;
    double margin = 0.2;
    double target = -2.904;
    if(abs(E-target)<=margin){
        B = true;
    }
    CHECK(B);
}



TEST(HFSolve){

    // Test HFSolve using the known value of the groundstate energy of Beryllium
    // as a control.



    double nProtons  = 4; //number of protons
    int nElectrons= 4; //number of electrons

    basis BS; //initialize basis object
    BS.init_HTO4(nProtons);


    hartreefocksolver object (BS,nElectrons,nProtons);  //initialize solver using 4 protons in the nucleus and 3 contracted orbitals
    double E = object.solve();                          //solve for the given basis


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
    AB.pp(primitiveA,primitiveB);
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

TEST(GaussianElectronElectron_test1){

    vec posA = {-0.5, 0, 0};
    vec posB = {-0.5, 0, 0};
    vec posC = {-0.5, 0, 0};
    vec posD = {-0.5, 0, 0};
    double a = 13.0077;
    double b = 13.0077;
    double c = 13.0077;
    double d = 13.0077;

    BoysFunction boys(2);

    Primitive primitiveA(1.0, 0, 0 ,0, a, posA);
    Primitive primitiveB(1.0, 0, 0 ,0, b, posB);
    Primitive primitiveC(1.0, 0, 0 ,0, c, posC);
    Primitive primitiveD(1.0, 0, 0 ,0, d, posD);

    integrator AB(primitiveA,primitiveB,boys);
    double particleParticleIntegral = AB.pp(primitiveC,primitiveD);
    cout << "----------------------------------------------------------------" << endl;
    cout << "-------------- TEST Particle-Particle Integral 1 ---------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << " pp = " << particleParticleIntegral << " and should be:" << endl;
    cout << " pp = 0.0071666040410096028615" << endl;
    cout << "________________________________________________________________" << endl;
    // regression test
    CHECK_CLOSE(0.0071666040410096028615, particleParticleIntegral, 1.0e-13);
}

TEST(GaussianElectronElectron_test2){
    vec posA = {0.5, 0, 0};
    vec posB = {-0.5, 0, 0};
    vec posC = {-0.5, 0, 0};
    vec posD = {0.5, 0, 0};
    double a = 13.0077;
    double b = 0.121949;
    double c = 0.444529;
    double d = 13.0077;

    BoysFunction boys(2);

    Primitive primitiveA(1.0, 0, 0 ,0, a,posA);
    Primitive primitiveB(1.0, 0, 0 ,0, b,posB);
    Primitive primitiveC(1.0, 0, 0 ,0, c,posC);
    Primitive primitiveD(1.0, 0, 0 ,0, d,posD);

    integrator AB(primitiveA,primitiveB,boys);
    double particleParticleIntegral = AB.pp(primitiveC,primitiveD);
    cout << "----------------------------------------------------------------" << endl;
    cout << "-------------- TEST Particle-Particle Integral 2 ---------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << " pp = " << particleParticleIntegral << " and should be:" << endl;
    cout << " pp = 0.022124581472837051566" << endl;
    cout << "________________________________________________________________" << endl;
    // regression test
    CHECK_CLOSE(0.022124581472837051566, particleParticleIntegral,1.0e-13);
}

TEST(GaussianElectronElectron_test3){
    vec posA = {0.5, 0, 0};
    vec posB = {-0.5, 0, 0};
    vec posC = {-0.5, 0, 0};
    vec posD = {0.5, 0, 0};
    double a = 13.0077;
    double b = 0.121949;
    double c = 0.444529;
    double d = 13.0077;

    BoysFunction boys(2);

    Primitive primitiveA(1.0, 0, 0 ,0, a,posA);
    Primitive primitiveB(1.0, 0, 1 ,0, b,posB);
    Primitive primitiveC(1.0, 0, 1 ,0, c,posC);
    Primitive primitiveD(1.0, 0, 0 ,0, d,posD);

    integrator AB(primitiveA,primitiveB,boys);
    double particleParticleIntegral = AB.pp(primitiveC,primitiveD);
    cout << "----------------------------------------------------------------" << endl;
    cout << "-------------- TEST Particle-Particle Integral 3 ---------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << " pp = " << particleParticleIntegral << " and should be:" << endl;
    cout << " pp = 0.0001385810300677682" << endl;
    cout << "________________________________________________________________" << endl;
    // regression test
    CHECK_CLOSE(0.0001385810300677682,particleParticleIntegral,1.0e-13);
}

TEST(GaussianElectronElectron_test4){
    vec posA = {1.2,2.3,3.4};
    vec posB = {-1.3,1.4,-2.4};
    vec posC = {2.3,0.9,3.2};
    vec posD = {5.0,1.9,1.2};
    double a = 0.2;
    double b = 0.3;
    double c = 0.4;
    double d = 0.1;

    BoysFunction boys(2);

    Primitive primitiveA(1.0, 0, 0 ,0, a,posA);
    Primitive primitiveB(1.0, 1, 0 ,0, b,posB);
    Primitive primitiveC(1.0, 0, 2 ,0, c,posC);
    Primitive primitiveD(1.0, 0, 0 ,1, d,posD);

    integrator AB(primitiveA,primitiveB,boys);
    double particleParticleIntegral = AB.pp(primitiveC,primitiveD);
    cout << "----------------------------------------------------------------" << endl;
    cout << "-------------- TEST Particle-Particle Integral 4 ---------------" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << setprecision(13) << " pp = " << particleParticleIntegral << " and should be:" << endl;
    cout << " pp = 0.2681206720738772" << endl;
    cout << "________________________________________________________________" << endl;
    // regression test
    CHECK_CLOSE(0.2681206720738772,particleParticleIntegral,1.0e-13);
}

int main() {

    return UnitTest::RunAllTests();
}
