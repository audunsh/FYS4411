#include "primitive.h"
#include <cmath>
#include <armadillo>

using namespace arma;

Primitive::Primitive(double weight, int xExponent, int yExponent, int zExponent, double Exponent, vec nucleusPosition) :
    m_weight(weight),
    m_xExponent(xExponent),
    m_yExponent(yExponent),
    m_zExponent(zExponent),
    m_Exponent(Exponent),
    m_nucleusPosition(nucleusPosition) {}


double Primitive::Exponent() const
{
    return m_Exponent;
}
int Primitive::zExponent () const
{
    return m_zExponent;
}
int Primitive::yExponent() const
{
    return m_yExponent;
}
int Primitive::xExponent()const
{
    return m_xExponent;
}
double Primitive::weight() const
{
    return m_weight;
}
const vec& Primitive::nucleusPosition() const
{
    return m_nucleusPosition;
}
