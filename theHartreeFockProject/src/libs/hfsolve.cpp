#include <hfsolve.h>
#include <lib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <basis.h>
#include <rhfsolve.h>
#include <uhfsolve.h>

using namespace arma;
using namespace std;

HFSolve::HFSolve(){}

HFSolve::HFSolve(basis BS){
    Bs = BS;
}

void HFSolve::solve_rhf(int N_electrons){
    //Do a restricted HF-calc and extract the variables needed by CCSolve
    rhfsolve object (Bs, N_electrons);

    energy = object.solve();
    C = object.C;
    epsilon = object.epsilon;
}

void HFSolve::solve_uhf(int N_electrons_up, int N_electrons_down){
    //Do a restricted HF-calc and extract the variables needed by CCSolve
    uhfsolve object (Bs, N_electrons_up, N_electrons_down);
    energy = object.solve();
    C = object.C;
    epsilon = object.epsilon;
}

void HFSolve::reset(){}



