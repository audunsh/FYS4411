#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <primitive.h>

class integrator
{
public:
    integrator();
    double overlap(primitive a, primitive b);
    void setupHermiteCoefficients();


};

#endif // INTEGRATOR_H
