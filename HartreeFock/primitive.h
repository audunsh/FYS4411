#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <vector>
#include <armadillo>

using namespace arma;

/* Primitive sets up a contracted basis function from
 * primitive gaussian functions.
 */

class Primitive {
public:
    explicit Primitive(double weight, int xExponent, int yExponent, int zExponent, double Exponent, vec nucleusPosition);

    double weight() const;
    int xExponent() const;
    int yExponent() const;
    int zExponent() const;
    double Exponent() const;
    const vec& nucleusPosition() const;

private:

    double m_weight;
    int m_xExponent;
    int m_yExponent;
    int m_zExponent;
    int m_Exponent;
    vec m_nucleusPosition;
};

#endif // PRIMITIVE_H
