#include "integrator.h"
#include <primitive.h>
#include <armadillo>

using namespace std;
using namespace arma;

integrator::integrator()
{
}

double integrator::overlap(primitive a, primitive b){
    return 0;
}

void integrator::setupHermiteCoefficients(){
}
