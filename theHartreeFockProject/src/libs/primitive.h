#ifndef PRIMITIVE H
#define PRIMITIVE H
#include <armadillo>
using namespace arma;

class Primitive {

public:
    explicit Primitive(double weight, int xExponent,
                       int yExponent, int zExponent, double exponent,
                       vec nucleusPosition);

    double exponent() const;
    int zExponent() const;
    int yExponent() const;
    int xExponent() const;
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

#endif // PRIMITIVE H
