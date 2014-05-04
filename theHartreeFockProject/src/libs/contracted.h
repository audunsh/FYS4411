#ifndef CONTRACTED_H
#define CONTRACTED_H
#include <armadillo>
#include <primitive.h>


using namespace std;
using namespace arma;

class contracted
{
public:
    contracted(int N,Primitive primitives[]);
    Primitive getPrimitive(int n);
    void setPrimitive(int n);
    void free();
private:
    Primitive basisFunction[];
    int Nprimitives;

};

#endif // CONTRACTED_H
