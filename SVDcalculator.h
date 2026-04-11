#ifndef VERIFYBAMID_SVDCALCULATOR_H
#define VERIFYBAMID_SVDCALCULATOR_H


#include <vector>
#include <map>
#include <unordered_set>
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

    /// Read a reference VCF, compute SVD, and write .UD, .mu, .bed, .V files.
    /// @param includeChr  if non-empty, only markers on these chromosomes are used
    /// @param skipMinSampleCountCheck  if true, < 1000 samples becomes a warning not an error
    void ProcessRefVCF(const std::string& VcfPath,
                       const std::unordered_set<std::string>& includeChr,
                       bool skipMinSampleCountCheck = false);

    /// Parse a VCF into a genotype matrix.
    /// @param includeChr  if non-empty, only markers on these chromosomes are used
    int ReadVcf(const std::string &VcfPath,
                std::vector<std::vector<char> >& genotype,
                int & nSamples, int& nMarkers,
                const std::unordered_set<std::string>& includeChr);

    std::vector<std::vector<double>> GetUDMatrix();
    std::vector<std::vector<double>> GetPCMatrix();
    std::vector<PCtype> GetMuArray();
    BED GetchooseBed();
    std::vector<region_t> GetBedVec();
    void WriteSVD(const std::string &VcfPath);
};


#endif //VERIFYBAMID_SVDCALCULATOR_H
