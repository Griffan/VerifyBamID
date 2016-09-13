#include <iostream>
#include "ContaminationEstimator.h"
#include "params.h"

using namespace std;

int main(int argc, char **argv) {
    cout << "Hello, World!" << endl;

    string UDPath("Empty"), MeanPath("Empty"), BedPath("Empty"), BamFileList("Empty"), BamFile("Empty"), RefPath(
            "Empty");
    string knownAF("Empty");
    bool fixPC(false), fixAlpha(false);
    int nfiles(0),seed(12345);
    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
                    LONG_PARAM_GROUP("Input/Output Files",
                                     "Input/Output files for the program[Complete Path Recommended]")
                    LONG_STRING_PARAM("UDPath", &UDPath,
                                      "[String] UD matrix file from SVD result of genotype matrix[Required]")
                    LONG_STRING_PARAM("MeanPath", &MeanPath,
                                      "[String] Mean matrix file of genotype matrix[Required]")
                    LONG_STRING_PARAM("BamFile", &BamFile,
                                      "[String] Bam or Cram file for the sample[Required]")
                    LONG_STRING_PARAM("BedPath", &BedPath,
                                      "[String] Bed file for markers used in this analysis,(chr\tpos-1\tpos\trefAllele\taltAllele)[Required]")
                    LONG_STRING_PARAM("Reference", &RefPath,
                                      "[String] reference file[Required]")
                    LONG_INT_PARAM("Seed",&seed,"[INT] Random number seed(default:12345)")
                    LONG_PARAM("fixPC", &fixPC, "[Bool] Fix PCs to estimate alpha[Optional]")
                    LONG_PARAM("fixAlpha", &fixAlpha, "[Bool] fixAlpha to estimate PC coordinates[Optional]")
                    LONG_STRING_PARAM("knownAF", &knownAF, "[String] known allele frequency file (chr\tpos\tfreq)[Optional]")



    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();
    if (UDPath == "Empty") {
        error("--UDPath is required");
        exit(EXIT_FAILURE);
    }
    if (BedPath == "Empty") {
        error("--BedPath is required");
        exit(EXIT_FAILURE);
    }
    if (MeanPath == "Empty") {
        error("--MeanPath is required");
        exit(EXIT_FAILURE);
    }
    if (RefPath == "Empty") {
        error("--Reference is required");
        exit(EXIT_FAILURE);
    }
    if (BamFile != "Empty") {
        nfiles = 1;
    }
//    else if (BamFileList !="Empty")
//    {
//        nfiles=2;
//    }
    else {
        error("--BamFile is required");
        exit(EXIT_FAILURE);
    }

    ContaminationEstimator Estimator(BamFile.c_str(), RefPath.c_str(), BedPath.c_str());
    Estimator.seed=seed;
    Estimator.ReadSVDMatrix(UDPath, MeanPath, BedPath);
    //std::cerr<<"NumMarker:"<<Estimator.NumMarker<<" and UD size:"<<Estimator.UD.size()<<std::endl;
    if(fixPC)
        Estimator.isPCFixed = true;
    else if (fixAlpha)
        Estimator.isAlphaFixed = true;
    if(knownAF!="Empty") {
        Estimator.isAFknown = true;
        Estimator.ReadAF(knownAF);
    }

    Estimator.OptimizeLLK();
    return 0;
}