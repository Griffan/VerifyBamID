//
// Created by Fan Zhang on 11/2/16.
//

#ifndef VERIFYBAMID_SVDCALCULATOR_H
#define VERIFYBAMID_SVDCALCULATOR_H


#include <vector>

class SVDcalculator {
private:
    int numIndividual;
    int numMarker;
#define PCtype double

    std::vector<std::vector<PCtype> > UD;
    std::vector<std::vector<PCtype> > PC;
    std::vector<PCtype> Mu;
public:
    SVDcalculator();
    ~SVDcalculator();
    int ReadVcf(const std::string& VcfPath);
    int Decompose();
    std::vector<std::vector<PCtype> >& GetUDMatrix();
    std::vector<std::vector<PCtype> >& GetPCMatrix();
    std::vector<PCtype>& GetMuArray();
};


#endif //VERIFYBAMID_SVDCALCULATOR_H
