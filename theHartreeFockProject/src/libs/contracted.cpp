#include "contracted.h"
#include <armadillo>
#include <primitive.h>

using namespace std;
using namespace arma;


contracted::contracted(int N)
{

    //lincomb.set_size(N);
}

void contracted::setPrimitive(int n, Primitive P){
    //lincomb(n) = P;
}

Primitive contracted::getPrimitive(int n){
    return Primitive(0,0,0,0,0,{0,0,0});//lincomb(n);
}
