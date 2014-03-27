#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <primitive.h>

class integrator
{
public:
    integrator();
    double overlapIntegral(Primitive a, Primitive b);
    void setupHermiteCoefficients();


};

#endif // INTEGRATOR_H
