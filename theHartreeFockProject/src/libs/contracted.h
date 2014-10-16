#ifndef CONTRACTED_H
#define CONTRACTED_H
#include <armadillo>
#include <primitive.h>


using namespace std;
using namespace arma;

class contracted
{
public:
    contracted();
    contracted(int N,Primitive primitives[]);
    Primitive getPrimitive(int n);
    void setPrimitive(int n);
    void appendPrimitive(Primitive P);
    void free();
private:
    Primitive basisFunction[];
    vector<Primitive> basisFs;
    int Nprimitives;

};

#endif // CONTRACTED_H
