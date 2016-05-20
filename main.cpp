#include <iostream>
#include "ContaminationEstimator.h"
#include "params.h"
using namespace std;

int main(int argc, char ** argv) {
    cout << "Hello, World!" << endl;

    string UDPath("Empty"),MeanPath("Empty"),BedPath("Empty"),BamFileList("Empty"),BamFile("Empty"),RefPath("Empty");
    int nfiles(0);
    paramList pl;
    BEGIN_LONG_PARAMS(longParameters) LONG_PARAM_GROUP("Input/Output Files", "Input/Output files for the program[Complete Path Recommended]")
    LONG_STRING_PARAM("UDPath", &UDPath, "[String] Pair end 1 fastq file[Leave empty if using fq_list or bam_in]")
//    LONG_STRING_PARAM("BamFileList", &BamFileList, "[String] Pair end 2 fastq file.[Leave empty if using single end]")
    LONG_STRING_PARAM("BamFile", &BamFile, "[String] Pair end 2 fastq file.[Leave empty if using single end]")
    LONG_STRING_PARAM("BedPath", &BedPath, "[String] Path of input fastq files, tab-delimited, one pair-end files per line(one file per line for single end)[Leave empty if using bam_in or fastq_1]")
    LONG_STRING_PARAM("MeanPath", &MeanPath, "[String] Input bam file path[Leave empty if using fq_list or fastq_1]")
    LONG_STRING_PARAM("Reference", &RefPath, "[String] Prefix of all the output files[Required]")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();
    if (UDPath == "Empty")
    {
        error("--UDPath is required");
        exit(EXIT_FAILURE);
    }
    if (BedPath == "Empty")
    {
        error("--BedPath is required");
        exit(EXIT_FAILURE);
    }
    if (MeanPath == "Empty")
    {
        error("--MeanPath is required");
        exit(EXIT_FAILURE);
    }
    if (RefPath == "Empty")
    {
        error("--Reference is required");
        exit(EXIT_FAILURE);
    }
    if (BamFile != "Empty")
    {
        nfiles=1;
    }
//    else if (BamFileList !="Empty")
//    {
//        nfiles=2;
//    }
    else
    {
        error("Either --BamFile is required");
        exit(EXIT_FAILURE);
    }

    ContaminationEstimator Estimator(BamFile.c_str(),RefPath.c_str(),BedPath.c_str(),nfiles);
    Estimator.ReadSVDMatrix(UDPath,MeanPath,BedPath);
    //std::cerr<<"NumMarker:"<<Estimator.NumMarker<<" and UD size:"<<Estimator.UD.size()<<std::endl;
    Estimator.OptimizeLLK();
    return 0;
}