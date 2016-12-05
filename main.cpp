/* The MIT License

   Copyright (c) 2009 Genome Research Ltd (GRL), 2010 Broad Institute

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include "ContaminationEstimator.h"
#include "params.h"
#include "SVDcalculator.h"

using namespace std;

int main(int argc, char **argv) {
    cout << "Hello, World!" << endl;

    string UDPath("Empty"), MeanPath("Empty"), BedPath("Empty"), BamFileList("Empty"), BamFile("Empty"), RefPath(
            "Empty"), OutputPrefix("result");
    string knownAF("Empty");
    string RefVCF("Empty");
    string fixPC("Empty");
    double fixAlpha(std::numeric_limits<double>::max());
    bool asHeter(false);
    int nfiles(0),seed(12345),nPC(2);
    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
                    LONG_PARAM_GROUP("Input/Output Files",
                                     "Input/Output files for the program[Complete Path Recommended]")
                    LONG_STRING_PARAM("BamFile", &BamFile,
                                      "[String] Bam or Cram file for the sample[Required]")
                    LONG_STRING_PARAM("BedPath", &BedPath,
                                      "[String] Bed file for markers used in this analysis,1 based pos(chr\tpos-1\tpos\trefAllele\taltAllele)[Required]")
                    LONG_STRING_PARAM("Reference", &RefPath,
                                      "[String] reference file[Required]")
                    LONG_STRING_PARAM("RefVCF", &RefVCF,
                                      "[String] VCF file from which to extract reference panel's genotype matrix")
                    LONG_STRING_PARAM("UDPath", &UDPath,
                                      "[String] UD matrix file from SVD result of genotype matrix")
                    LONG_STRING_PARAM("MeanPath", &MeanPath,
                                      "[String] Mean matrix file of genotype matrix")
                    LONG_STRING_PARAM("Output", &OutputPrefix,
                                      "[String] OutputPrefix[optional]")
                    LONG_INT_PARAM("numPC", &nPC,
                                   "[Int] number of PCs used to infer Allele Frequency[optional]")
                    LONG_STRING_PARAM("fixPC", &fixPC, "[String] Input fixed PCs to estimate Alpha[format:PC1|PC2|PC3...]")
                    LONG_DOUBLE_PARAM("fixAlpha", &fixAlpha, "[Double] Input fixed Alpha to estimate PC coordinates")
//                    LONG_PARAM("fixPC", &fixPC, "[Bool] Fix PCs to estimate localAlpha[default:false]")
//                    LONG_PARAM("fixAlpha", &fixAlpha, "[Bool] fixAlpha to estimate PC coordinates[default:false]")
                    LONG_PARAM("asHeter", &asHeter, "[Bool] Infer contamination level as if target sample and contamination source are from the different population[default:false]")
                    LONG_STRING_PARAM("knownAF", &knownAF, "[String] known allele frequency file (chr\tpos\tfreq)[Optional]")
                    LONG_INT_PARAM("Seed",&seed,"[INT] Random number seed[default:12345]")

    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if (RefVCF == "Empty") {
        if (UDPath == "Empty") {
            error("--UDPath is required when --RefVCF is absent");
            exit(EXIT_FAILURE);
        }
        if (MeanPath == "Empty") {
            error("--MeanPath is required when --RefVCF is absent");
            exit(EXIT_FAILURE);
        }
        if (BedPath == "Empty") {
            error("--BedPath is required when --RefVCF is absent");
            exit(EXIT_FAILURE);
        }
    } else//SVD on the fly
    {
        notice("Specified reference panel VCF file, doing SVD on the fly...");
        notice("This procedure will generate SVD matrices as [RefVCF path].UD and [RefVCF path].mu");
        notice("You may specify --UDPath [RefVCF path].UD and --MeanPath [RefVCF path].mu in future use");
        SVDcalculator calculator;
        calculator.ProcessRefVCF(RefVCF);
        UDPath = RefVCF+".UD";
        MeanPath = RefVCF+".mu";
        BedPath = RefVCF+".bed";
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

    ContaminationEstimator Estimator(nPC, BamFile.c_str(), RefPath.c_str(), BedPath.c_str());
    Estimator.seed = seed;
    Estimator.isHeter = asHeter;
    Estimator.ReadSVDMatrix(UDPath, MeanPath);
    //std::cerr<<"NumMarker:"<<Estimator.NumMarker<<" and UD size:"<<Estimator.UD.size()<<std::endl;
    if(fixPC!="Empty") {// parse --fixPC
        notice("you specified --fixPC, this will overide dynamic estimation of PCs");
        notice("parsing the PCs");
        stringstream ss(fixPC);
        string token;
        std::vector<PCtype> tmpPC;
        while(std::getline(ss, token, '|')) {
            tmpPC.push_back(atof(token.c_str()));
        }
        if(tmpPC.size() > nPC)
        {
            warning("parameter --fixPC provided larger dimension than parameter --numPC(default value 2) and hence will be truncated");
            for (int i = 0; i < nPC; ++i) {
                Estimator.PC[0][i]=tmpPC[i];
            }
        }
        Estimator.isPCFixed = true;
    }
    else if(abs(fixAlpha - std::numeric_limits<double>::max())>std::numeric_limits<double>::epsilon()) {
        notice("you specified --fixAlpha, this will overide dynamic estimation of alpha");
        Estimator.alpha = fixAlpha;
        Estimator.isAlphaFixed = true;
    }
    if(knownAF!="Empty") {
        Estimator.isAFknown = true;
        Estimator.ReadAF(knownAF);
    }

    Estimator.OptimizeLLK();

    const char* headers[] = {"#SEQ_ID","RG","CHIP_ID","#SNPS","#READS","AVG_DP","FREEMIX","FREELK1","FREELK0","FREE_RH","FREE_RA","CHIPMIX","CHIPLK1","CHIPLK0","CHIP_RH","CHIP_RA","DPREF","RDPHET","RDPALT"};
    std::ofstream fout(OutputPrefix+".selfSM");
    fout<<headers<<std::endl;


    return 0;
}
