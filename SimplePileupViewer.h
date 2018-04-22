#ifndef CONTAMINATIONFINDER_SIMPLEPILEUP_H
#define CONTAMINATIONFINDER_SIMPLEPILEUP_H
#include "sam_opts.h"
#include "htslib/faidx.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>

typedef struct {
    int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag;
    int rflag_require, rflag_filter;
    int openQ, extQ, tandemQ, min_support; // for indels
    double min_frac; // for indels
    char *reg, *pl_list, *fai_fname, *output_fname;
    faidx_t *fai;
    void *bed, *rghash;
    int argc;
    char **argv;
    sam_global_args ga;
} mplp_conf_t;


typedef std::vector<std::vector<char> >BaseInfo;
typedef std::vector<std::vector<char> >QualInfo;

class ContaminationEstimator;
class SVDcalculator;
struct region_t{
    std::string chr;
    int beg;//0 based
    int end;
    region_t(std::string chr0, int beg0, int end0)
    {
        chr=chr0;
        beg=beg0;
        end=end0;
    }
    region_t()
    {
        chr="";
        beg=0;
        end=0;
    }
} ;

typedef std::unordered_map<std::string, std::unordered_map<int, std::pair<char, char> > > BED;

class SimplePileupViewer {

public:

    int regIndex;
    std::vector<region_t>* bedVec;
    BED bedTable;

    mplp_conf_t mplp;

    BaseInfo baseInfo;
    QualInfo qualInfo;

    std::string SEQ_SM = "DefaultSampleName";
    int numBases = 0;
    int effectiveNumSite = 0;
    double avgDepth = 0;
    double sdDepth = 0;

    double firstQT;//1st quantile
    double thirdQT;//2st quantile

    std::unordered_map<std::string,std::unordered_map<int32_t,int32_t> > posIndex;

    SimplePileupViewer();

    SimplePileupViewer(std::vector<region_t>* A,const char *bamFile, const char* faiFile, const char* bedFile,int nfiles=1);

    int SIMPLEmpileup(mplp_conf_t *conf, int n, char **fn);

    SimplePileupViewer(const BED& BedFromEstimator, const std::string &pileupFile);

    int ReadPileup(const std::string & filePath);

    region_t GetNextRegion(bool& BedEOF)
    {
        BedEOF=false;
        if(regIndex<bedVec->size()) {
            if(regIndex==(bedVec->size()-1)) BedEOF=true;
            return bedVec->at(regIndex++);
        }
        else {
            std::cerr<<"Region number out of bound!\n"<<std::endl;
            exit(EXIT_FAILURE);
        }
    }

    int GetNumMarker()
    {
        return effectiveNumSite;
    }

    inline std::vector<char>& GetBaseInfoAt(std::string& chr, int32_t pos)
    {
        return baseInfo[posIndex[chr][pos]];
    }
    inline std::vector<char>& GetQualInfoAt(std::string& chr, int32_t pos)
    {
        return qualInfo[posIndex[chr][pos]];
    }

    virtual ~SimplePileupViewer();
};


#endif //CONTAMINATIONFINDER_SIMPLEPILEUP_H
