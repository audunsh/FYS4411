#include "contracted.h"
#include <armadillo>
#include <primitive.h>


using namespace std;
using namespace arma;


contracted::contracted(int N, Primitive primitives[]){
    basisFunction[N];
    for(int i=0; i<N; i++){
        basisFunction[i] = primitives[i];
    }
}

void contracted::setPrimitive(int n){
}

Primitive contracted::getPrimitive(int n){
    return basisFunction[n];
}

void contracted::free(){
    delete[] basisFunction;
}
