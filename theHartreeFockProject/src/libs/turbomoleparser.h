#ifndef TURBOMOLEPARSER_H
#define TURBOMOLEPARSER_H

#include <string>
#include <unordered_map>
#include <vector>
#include <basis.h>

using namespace std;

class TurboMoleParser
{
public:
    TurboMoleParser();

    bool load(string fileName);

    //const vector<GaussianContractedOrbital> &contractedBasisFunctions() const;

    //HF::AtomType atomType() const;

    double normalizationFactor(double exp, urowvec pows);
    int factorial(int n);

private:
    //vector<GaussianContractedOrbital> m_contractedBasisFunctions;
    //HF::AtomOrbitalType m_currentOrbitalType = HF::sOrbitalType;
    //vector<GaussianPrimitiveOrbital> m_collectedPrimitiveBasisFunctions;
    //void mergePrimitivesIntoContracted();
    //HF::AtomType m_atomType;
};

#endif // TURBOMOLEPARSER_H
