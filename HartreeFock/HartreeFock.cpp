/* Computational Physics II FYS4411 - the Hartree-Fock project - Goeran & Audun - Masterstudents at CompPhys UIO
 *
 * The Hartree-Fock method:
 * is a variational method in which the wave functions of a many-electron-system take the form of an
 * antisymmetised product of one-electron wave functions. This restriction leads to an effectice
 * Schroredinger equation for the individual one-electron wave functions (orbitals) with a potential
 * determined by the orbitals occupied by the other electrons. The coupling between the orbitals via
 * the potentials causes the resulting equations to become non-linear in the orbitals, and the
 * solution must be found iteratively in a self-consistency procedure. The Hartree-Fock method is a
 * procedure close to a statistical mechanics approach.
 */

using namespace std;
//using namespace arma;

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
//#include <lib.h>

int main(){
    cout << "Hello to everyone looking at the screen!" << endl;
    return 0;
} // End: main()
