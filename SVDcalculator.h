#ifndef VERIFYBAMID_SVDCALCULATOR_H
#define VERIFYBAMID_SVDCALCULATOR_H


#include <vector>
#include <map>
#include "SimplePileupViewer.h"

class SVDcalculator {
private:
    int numIndividual;
    int numMarker;
#define PCtype double

    std::vector<std::vector<PCtype> > UD;
    std::vector<std::vector<PCtype> > PC;
    std::vector<std::string> Samples;
    std::vector<PCtype> Mu;
    BED chooseBed;
    std::vector<region_t> BedVec;
public:
    SVDcalculator();
    ~SVDcalculator();
    void ProcessRefVCF(const std::string& VcfPath);
    int ReadVcf(const std::string &VcfPath,
                std::vector<std::vector<char> >& genotype,
                int & nSamples, int& nMarkers);
//    int Decompose();
    std::vector<std::vector<double>> GetUDMatrix();//return the matrix
    std::vector<std::vector<double>> GetPCMatrix();
    std::vector<PCtype> GetMuArray();
    BED GetchooseBed();
    std::vector<region_t> GetBedVec();
    void WriteSVD(const std::string &VcfPath);
};


#endif //VERIFYBAMID_SVDCALCULATOR_H
