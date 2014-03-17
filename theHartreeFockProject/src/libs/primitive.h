#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <armadillo>

using namespace arma;

class primitive
{
public:
    explicit primitive(double weight, int xExponent, int yExponent, int zExponent, double exponent, vec nucleusPosition);
    double exponent() const;
    int xExponent() const;
    int yExponent() const;
    int zExponent() const;
    double weight() const;
    const vec& nucleusPosition() const;
private:
    double m_weight;
    int m_xExponent;
    int m_yExponent;
    int m_zExponent;
    double m_exponent;
    vec m_nucleusPosition;
};

#endif // PRIMITIVE_H
