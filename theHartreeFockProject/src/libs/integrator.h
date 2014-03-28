#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <primitive.h>

class integrator
{
public:
    integrator();
    double overlapIntegral(Primitive a, Primitive b);
    void setupHermiteCoefficients();

private:
    double pi = 4*atan(1);

};

#endif // INTEGRATOR_H
