#include "contracted.h"
#include <armadillo>
#include <primitive.h>


using namespace std;
using namespace arma;


contracted::contracted(int N, Primitive primitives[]){
    Nprimitives = N;
    //basisFunction[N];
    //vector basisFs.resize(N);
    //basisFs.resize(N);
    for(int i=0; i<N; i++){
        //basisFunction[i] = primitives[i];
        basisFs.push_back(primitives[i]);
    }
}

void contracted::setPrimitive(int n){
}

Primitive contracted::getPrimitive(int n){
    //return basisFunction[n];
    return basisFs[n];
}

void contracted::free(){
    delete[] basisFunction;
}
