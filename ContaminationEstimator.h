/*The MIT License (MIT)

Copyright (c) 2017 Fan Zhang, Hyun Min Kang

        Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */
/* Contact: Fan Zhang <fanzhang@umich.edu> */
#ifndef CONTAMINATIONESTIMATOR_H_
#define CONTAMINATIONESTIMATOR_H_

#include <string>
#include <unordered_map>
//#include <tkDecls.h>
#include "MathVector.h"
#include "MathGenMin.h"
#include "SimplePileupViewer.h"
#include <limits>

#ifdef _OPENMP
#include "omp.h"
#endif

class ContaminationEstimator {
public:
    bool isPCFixed;
    bool isAlphaFixed;
    bool isAFknown;
    bool isHeter;
    bool isPileupInput;
    bool isSanityCheckDisabled;
    bool verbose;
    int numPC;
    int numThread;
    int seed;
    double epsilon;
#define PCtype double

//#define PHRED(x)    pow(10.0,x/-10.0)
    static double Phred(double x) {
        return pow(10.0, x / -10.0);
    }

    class FullLLKFunc : public VectorFunc {
    public:
        double min_af;
        double max_af;
        double llk1;
        double llk0;
        ContaminationEstimator *ptr;
        std::vector<double> fixPC;
        std::vector<double> fixPC2;
        double fixAlpha;
        std::vector<double> globalPC;//best result holder
        std::vector<double> globalPC2;//best result holder
        double globalAlpha;//best result holder
        const char *Base;

        FullLLKFunc() {
            FullLLKFunc::Base = "actg";
            min_af = 0.00005;
            max_af = 0.99995;
            llk1 = 0;
            ptr = nullptr;
            fixAlpha = 0;
            std::cerr << "Initialize from FullLLKFunc()" << std::endl;

        }

        FullLLKFunc(int dim, ContaminationEstimator *contPtr) : fixPC(dim, 0.), fixPC2(dim, 0.), globalPC(fixPC),
                                                                globalPC2(fixPC2) {
            FullLLKFunc::Base = "actg";
            min_af = 0.00005;
            max_af = 0.99995;
            llk1 = 0.;
            ptr = contPtr;
            fixAlpha = 0.;
            globalAlpha = 0.;
            std::cerr << "Initialize from FullLLKFunc(int dim, ContaminationEstimator* contPtr)" << std::endl;
        }

        ~FullLLKFunc() {};

        inline static double InvLogit(double &x) {
            double e = exp(x);
            return e / (1. + e);
        };

        inline static double Logit(double &x) {

            return log(x / (1. - x));
        };

        inline int Normalize(std::vector<double> &tPC) {
            for (int i = 0; i < tPC.size(); ++i) {
                tPC[i] = (tPC[i] - ptr->muv[i]) / ptr->sdv[i];
            }
            return 0;
        };

        inline int InvNormalize(std::vector<double> &tPC) {
            for (int i = 0; i < tPC.size(); ++i) {
                tPC[i] = tPC[i] * ptr->sdv[i] + ptr->muv[i];
            }
            return 0;
        };

        inline char findAlt(std::vector<char> &tmpBase) {
            int a[4];
            int maxIndex(-1);
            for (int i = 0; i < tmpBase.size(); ++i) {
                if (tmpBase[i] == '.' || tmpBase[i] == ',') continue;
                if (tmpBase[i] == 'A' || tmpBase[i] == 'a') a[0]++;
                else if (tmpBase[i] == 'C' || tmpBase[i] == 'c') a[1]++;
                else if (tmpBase[i] == 'T' || tmpBase[i] == 't') a[2]++;
                else if (tmpBase[i] == 'G' || tmpBase[i] == 'g') a[3]++;
                maxIndex = 0;
            }
            if (maxIndex == -1) return 0;

            for (int j = 0; j < 4; ++j) {
                if (a[j] > a[maxIndex]) maxIndex = j;
            }
            return Base[maxIndex];
        }

        inline double getConditionalBaseLK(char base, int genotype, char altBase, bool is_error) {
            if (!is_error) {
                if (genotype == 0) {
                    if (base == '.' || base == ',') {
                        return 1;
                    } else
                        return 0;
                } else if (genotype == 1) {
                    if (base == '.' || base == ',') {
                        return 0.5;
                    } else if (toupper(base) == toupper(altBase)) {
                        return 0.5;
                    } else
                        return 0;
                } else if (genotype == 2) {
                    if (toupper(base) == toupper(altBase)) {
                        return 1;
                    } else
                        return 0;
                } else {
                    std::cerr << "genotype error!" << std::endl;
                    exit(EXIT_FAILURE);
                }
            } else {
                if (genotype == 0) {
                    if (base == '.' || base == ',') {
                        return 0;
                    } else if (toupper(base) == toupper(altBase)) {
                        return 1. / 3.;
                    } else
                        return 2. / 3.;
                } else if (genotype == 1) {
                    if (base == '.' || base == ',') {
                        return 1. / 6.;
                    } else if (toupper(base) == toupper(altBase)) {
                        return 1. / 6.;
                    } else
                        return 2. / 3.;
                } else if (genotype == 2) {
                    if (base == '.' || base == ',') {
                        return 1. / 3.;
                    }
                    if (toupper(base) == toupper(altBase)) {
                        return 0;
                    } else
                        return 2. / 3.;
                } else {
                    std::cerr << "genotype error!" << std::endl;
                    exit(EXIT_FAILURE);
                }

            }


        }

        void InitialGF(double AF, double *GF) const {
            if (AF < min_af) AF = min_af;
            if (AF > max_af) AF = max_af;
            GF[0] = (1 - AF) * (1 - AF);
            GF[1] = 2 * (AF) * (1 - AF);
            GF[2] = AF * AF;
        }

        inline double
        ComputeMixLLKs(const std::vector<double> &tPC1, const std::vector<double> &tPC2, const double alpha) {

            double sumLLK(0);
#ifdef _OPENMP
            omp_set_num_threads(ptr->numThread);
#pragma omp parallel for reduction (+:sumLLK)
#endif
            for (size_t i = 0; i < ptr->NumMarker; ++i) {

                std::string chr = ptr->PosVec[i].first;
                int pos = ptr->PosVec[i].second;
                if (ptr->viewer.posIndex.find(chr) == ptr->viewer.posIndex.end()) {
                    continue;
                } else if (ptr->viewer.posIndex[chr].find(pos) == ptr->viewer.posIndex[chr].end()) {
                    continue;
                }

                std::vector<char> tmpBase = ptr->viewer.GetBaseInfoAt(chr, pos);
                std::vector<char> tmpQual = ptr->viewer.GetQualInfoAt(chr, pos);

                if (tmpBase.size() == 0) continue;

                if (not ptr->isSanityCheckDisabled and
                    (tmpBase.size() < (ptr->viewer.avgDepth - 3 * ptr->viewer.sdDepth) or
                     tmpBase.size() > (ptr->viewer.avgDepth + 3 * ptr->viewer.sdDepth)))
                    continue;

                if (ptr->isAFknown) {
                    ptr->AFs[i] = ptr->AF2s[i] = ptr->knownAF[chr][pos];
                } else {
                    ptr->AFs[i] = 0.;
                    for (int k = 0; k < tPC1.size(); ++k) {
                        ptr->AFs[i] += ptr->UD[i][k] * tPC1[k];
                    }
                    ptr->AFs[i] += ptr->means[i];
                    ptr->AFs[i] /= 2.0;

                    ptr->AF2s[i] = 0.;
                    for (int k = 0; k < tPC2.size(); ++k) {
                        ptr->AF2s[i] += ptr->UD[i][k] * tPC2[k];
                    }
                    ptr->AF2s[i] += ptr->means[i];
                    ptr->AF2s[i] /= 2.0;
                }

                double markerLK(0);
                double GF[3];
                double GF2[3];

                InitialGF(ptr->AFs[i], GF);
                InitialGF(ptr->AF2s[i], GF2);

                char altBase = ptr->ChooseBed[chr][pos].second;

                for (int geno1 = 0; geno1 < 3; ++geno1)
                    for (int geno2 = 0; geno2 < 3; ++geno2) {
                        double baseLK(0);
                        for (int j = 0; j < tmpBase.size(); ++j) {
                            baseLK += log((alpha * getConditionalBaseLK(tmpBase[j], geno1, altBase, 1) +
                                           (1. - alpha) * getConditionalBaseLK(tmpBase[j], geno2, altBase, 1)) *
                                          Phred(tmpQual[j] - 33)
                                          + (alpha * getConditionalBaseLK(tmpBase[j], geno1, altBase, 0) +
                                             (1. - alpha) * getConditionalBaseLK(tmpBase[j], geno2, altBase, 0)) *
                                            (1 - Phred(tmpQual[j] - 33)));
//                            std::cerr <<i<<"th marker\t"<<tmpBase[j]<<"\t"<<tmpQual[j]<<"\t"<<altBase<<"\tlocalAlpha:"<<localAlpha<<"\tgeno1:"<<geno1<<"\tgeno2:"<<geno2
//                            <<"\tgetConditionalBaseLK1:"<<getConditionalBaseLK(tmpBase[j], geno1, altBase, 1)<<"\t"<< getConditionalBaseLK(tmpBase[j], geno2, altBase, 1)<<"\tPhred:"<<Phred(tmpQual[j] - 33)
//                            <<"\tgetConditionalBaseLK0:"<<getConditionalBaseLK(tmpBase[j], geno1, altBase, 0)<<"\t"<<getConditionalBaseLK(tmpBase[j], geno2, altBase, 0)<< std::endl;
                        }
                        markerLK += exp(baseLK) * GF[geno1] * GF2[geno2];
                    }
                if (markerLK > 0)
                    sumLLK += log(markerLK);
            }
            return sumLLK;
        }

        int Initialize() {
            globalPC = fixPC = globalPC2 = fixPC2 = ptr->PC[1];//only intended smaple has pre defined PCs
            globalAlpha = fixAlpha = ptr->alpha;
            llk1 = (0 - ComputeMixLLKs(fixPC, fixPC2, fixAlpha));

            for (int k = 0; k < ptr->numPC; ++k) {
                //ptr->PC[0][k] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
                ptr->PC[0][k] = 0.01;
            }
            for (int k = 0; k < ptr->numPC; ++k) {
                //ptr->PC[1][k] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
                ptr->PC[1][k] = 0.01;
            }
            //ptr->alpha = fabs(static_cast <double> (rand()) / static_cast <double> (RAND_MAX));
            ptr->alpha = 0.03;
            return 0;
        }

        int CalculateLLK0() {
            llk0 = (0 - ComputeMixLLKs(globalPC, globalPC, 0));
            return 0;
        }

        virtual double Evaluate(Vector &v) {
            double smLLK = 0;
            if (!ptr->isHeter) {
                if (ptr->isPCFixed) {
                    double tmpAlpha = InvLogit(v[0]);
                    smLLK = 0 - ComputeMixLLKs(fixPC, fixPC2, tmpAlpha);
                    if (smLLK < llk1) {
                        llk1 = smLLK;
                        globalAlpha = tmpAlpha;
                    }
                } else if (ptr->isAlphaFixed) {
                    std::vector<double> tmpPC(ptr->numPC, 0.);
                    for (int i = 0; i < ptr->numPC; ++i) {
                        tmpPC[i] = v[i];
                    }

                    smLLK = 0 - ComputeMixLLKs(tmpPC, tmpPC, fixAlpha);
                    if (smLLK < llk1) {
                        llk1 = smLLK;
                        globalPC = tmpPC;
                        globalPC2 = tmpPC;
                    }
                } else {
                    std::vector<double> tmpPC(ptr->numPC, 0.);
                    for (int i = 0; i < ptr->numPC; ++i) {
                        tmpPC[i] = v[i];
                    }
                    double tmpAlpha = InvLogit(v[ptr->numPC]);
                    smLLK = 0 - ComputeMixLLKs(tmpPC, tmpPC, tmpAlpha);
                    if (smLLK < llk1) {
                        llk1 = smLLK;
                        globalPC = tmpPC;
                        globalPC2 = tmpPC;
                        globalAlpha = tmpAlpha;
                    }
                }
            } else//contamination source from different population
            {
                if (ptr->isPCFixed) {//only fixed for intended sample
                    std::vector<double> tmpPC(ptr->numPC, 0.);
                    for (int i = 0; i < ptr->numPC; ++i) {
                        tmpPC[i] = v[i];
                    }
                    double tmpAlpha = InvLogit(v[ptr->numPC]);
                    smLLK = 0 - ComputeMixLLKs(tmpPC, fixPC2, tmpAlpha);

                    if (smLLK < llk1) {
                        llk1 = smLLK;
                        globalPC = tmpPC;
                        globalAlpha = tmpAlpha;
                    }
                } else if (ptr->isAlphaFixed) {
                    std::vector<double> tmpPC(ptr->numPC, 0.);
                    std::vector<double> tmpPC2(ptr->numPC, 0.);

                    for (int k = 0; k < v.Length(); ++k) {
                        if (k < ptr->numPC)
                            tmpPC[k] = v[k];
                        else if (k < ptr->numPC * 2)
                            tmpPC2[k - (ptr->numPC)] = v[k];
                        else {
                            error("Simplex Vector dimension error!");
                            exit(EXIT_FAILURE);
                        }
                    }
                    smLLK = 0 - ComputeMixLLKs(tmpPC, tmpPC2, fixAlpha);
                    if (smLLK < llk1) {
                        llk1 = smLLK;
                        globalPC = tmpPC;
                        globalPC2 = tmpPC2;
                    }
                } else {
                    std::vector<double> tmpPC(ptr->numPC, 0.);
                    std::vector<double> tmpPC2(ptr->numPC, 0.);
                    double tmpAlpha(0.);
                    for (int k = 0; k < v.Length(); ++k) {
                        if (k < ptr->numPC)
                            tmpPC[k] = v[k];
                        else if (k < ptr->numPC * 2)
                            tmpPC2[k - (ptr->numPC)] = v[k];
                        else if (k == ptr->numPC * 2)
                            tmpAlpha = InvLogit(v[k]);
                        else {
                            error("Simplex Vector dimension error!");
                            exit(EXIT_FAILURE);
                        }
                    }
                    smLLK = (0 - ComputeMixLLKs(tmpPC, tmpPC2, tmpAlpha));
                    if (smLLK < llk1) {
                        llk1 = smLLK;
                        globalPC = tmpPC;
                        globalPC2 = tmpPC2;
                        globalAlpha = tmpAlpha;
                    }
                }
            }
            if (ptr->verbose)
                std::cerr << "globalPC:" << globalPC[0] << "\tglobalPC:" << globalPC[1]
                          << "\tglobalPC2:" << globalPC2[0] << "\tglobalPC2:" << globalPC2[1]
                          << "\tglobalAlpha:" << globalAlpha << "\tllk:" << llk1 << std::endl;
            return smLLK;
        }
    };

    SimplePileupViewer viewer;
    uint32_t NumMarker;
    FullLLKFunc fn;

    std::unordered_map<std::string, std::unordered_map<uint32_t, double> > knownAF;

    double alpha;//input alpha
    std::vector<std::vector<PCtype> > UD;//input UD
    std::vector<std::vector<PCtype> > PC;//input PC
    std::vector<PCtype> means;
    ////
    std::vector<PCtype> muv;
    std::vector<PCtype> sdv;
    ////
    std::vector<double> AFs;
    std::vector<double> AF2s;

    typedef std::unordered_map<std::string, std::unordered_map<int, std::pair<char, char> > > BED;
    BED ChooseBed;//pos is 1-based
    std::vector<region_t> BedVec;//serialized BED info, convenient for bam reading
    std::vector<std::pair<std::string, int> > PosVec;

    ContaminationEstimator();

    ContaminationEstimator(int nPC, const char *bedFile, int nThread, double ep);

    /*Initialize from existed UD*/
    /*This assumes the markers are the same as the selected vcf*/
    /*ContaminationEstimator(const std::string &UDpath, const std::string &PCpath, const std::string &Mean,
                           const std::string &pileup, const std::string &GLpath, const std::string &Bed);
    */
    int ReadMatrixUD(const std::string &path);

    int ReadMatrixPC(const std::string &path);

    /*Intersect marker sites*/
    /*
    int ReadMatrixGL(const std::string &path);
    */
    int ReadChooseBed(const std::string &path);

    int ReadMean(const std::string &path);

    int ReadAF(const std::string &path);

    int ReadBam(const char *bamFile, const char *faiFile, const char *bedFile);

    int ReadPileup(const std::string &pileupFile);

    bool IsSanityCheckOK();
    /*
    int CheckMarkerSetConsistency();

    int FormatMarkerIntersection();
    */
    /*Optimize*/
    int OptimizeLLK(const std::string &OutputPrefix);

    ~ContaminationEstimator();

    /*
    int RunFromVCF(const std::string VcfSiteAFFile, const std::string CurrentMPU, const std::string ReadGroup,
                   const std::string Prefix);

    int RunFromSVDMatrix(const std::string UDpath, const std::string PCpath, const std::string Mean,
                         const std::string &MPUpath, const std::string &Bed, const std::string &Prefix,
                         const std::string &ReadGroup);
    */
    int ReadSVDMatrix(const std::string &UDpath, const std::string &PCpath, const std::string &Mean);

    /*
    int FromBamtoPileup();
     */
    bool OptimizeHomoFixedPC(AmoebaMinimizer &myMinimizer);

    bool OptimizeHomoFixedAlpha(AmoebaMinimizer &myMinimizer);

    bool OptimizeHomo(AmoebaMinimizer &myMinimizer);

    bool OptimizeHeterFixedPC(AmoebaMinimizer &myMinimizer);

    bool OptimizeHeterFixedAlpha(AmoebaMinimizer &myMinimizer);

    bool OptimizeHeter(AmoebaMinimizer &myMinimizer);
};

#endif /* CONTAMINATIONESTIMATOR_H_ */
