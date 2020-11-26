#ifndef CONTAMINATIONFINDER_SIMPLEPILEUP_H
#define CONTAMINATIONFINDER_SIMPLEPILEUP_H
#include "htslib/faidx.h"
#include "sam_opts.h"
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#define B2B_FMT_DP (1 << 0)
#define B2B_FMT_SP (1 << 1)
#define B2B_FMT_DV (1 << 2)

#define MPLP_BCF 1
#define MPLP_VCF (1 << 1)
#define MPLP_NO_COMP (1 << 2)
#define MPLP_NO_ORPHAN (1 << 3)
#define MPLP_REALN (1 << 4)
#define MPLP_NO_INDEL (1 << 5)
#define MPLP_REDO_BAQ (1 << 6)
#define MPLP_ILLUMINA13 (1 << 7)
#define MPLP_IGNORE_RG (1 << 8)
#define MPLP_PER_SAMPLE (1 << 9)
#define MPLP_SMART_OVERLAPS (1 << 10)

#define MPLP_PRINT_MAPQ_CHAR (1 << 11)
#define MPLP_PRINT_QPOS (1 << 12)
#define MPLP_PRINT_QNAME (1 << 13)
#define MPLP_PRINT_FLAG (1 << 14)
#define MPLP_PRINT_RNAME (1 << 15)
#define MPLP_PRINT_POS (1 << 16)
#define MPLP_PRINT_MAPQ (1 << 17)
#define MPLP_PRINT_CIGAR (1 << 18)
#define MPLP_PRINT_RNEXT (1 << 19)
#define MPLP_PRINT_PNEXT (1 << 20)
#define MPLP_PRINT_TLEN (1 << 21)
#define MPLP_PRINT_SEQ (1 << 22)
#define MPLP_PRINT_QUAL (1 << 23)

#define MPLP_MAX_DEPTH 8000
#define MPLP_MAX_INDEL_DEPTH 250

typedef struct {
  int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag,
      all, rev_del;
  int rflag_require, rflag_filter;
  int openQ, extQ, tandemQ, min_support; // for indels
  double min_frac;                       // for indels
  char *reg, *pl_list, *fai_fname, *output_fname;
  faidx_t *fai;
  void *bed, *rghash, *auxlist;
  int argc;
  char **argv;
  char sep, empty;
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

    SimplePileupViewer(std::vector<region_t> *A, const char *bamFile,
                       const char *faiFile, const char *bedFile,
                       mplp_conf_t *mplpPtr, int nfiles = 1);

    int SimplePileup(mplp_conf_t *conf, int n, char **fn);

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
