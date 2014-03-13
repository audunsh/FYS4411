#include "primitive.h"

primitive::primitive(double weight, int xExponent, int yExponent, int zExponent, double exponent, vec nucleusPosition) :
    m_weight(weight),
    m_xExponent(xExponent),
    m_yExponent(yExponent),
    m_zExponent(zExponent),
    m_exponent(exponent),
    m_nucleusPosition(nucleusPosition)

{
}

double primitive::exponent() const{
    return m_exponent;
}

int primitive::xExponent() const{
    return m_xExponent;
}
int primitive::yExponent() const{
    return m_yExponent;
}
int primitive::zExponent() const{
    return m_zExponent;
}

double primitive::weight() const{
    return m_weight;
}
const vec& primitive::nucleusPosition() const{
    return m_nucleusPosition;
}
