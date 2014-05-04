#ifndef CONTRACTED_H
#define CONTRACTED_H
#include <armadillo>
#include <primitive.h>

using namespace std;
using namespace arma;

class contracted
{
public:
    contracted(int N);
    Primitive getPrimitive(int n);
    void setPrimitive(int n, Primitive P);
private:
    Primitive lincomb [];

};

#endif // CONTRACTED_H
