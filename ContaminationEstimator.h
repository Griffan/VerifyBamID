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

#include <array>
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

    // Precomputed Phred error probabilities for quality scores 0-93 (after ASCII-33 offset).
    // Avoids repeated pow() calls in the hot loop.  C++11 guarantees that
    // initialization of a local static is thread-safe, so using an
    // immediately-invoked lambda ensures the table is built exactly once even
    // under concurrent first callers.
    static const double* getPhredTable() {
        static const std::array<double, 94> table = []() {
            std::array<double, 94> t{};
            for (int i = 0; i < 94; ++i) {
                t[i] = pow(10.0, i / -10.0);
            }
            return t;
        }();
        return table.data();
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
        const double *phredTable;

        FullLLKFunc() {
            FullLLKFunc::Base = "actg";
            min_af = 0.00005;
            max_af = 0.99995;
            llk1 = 0;
            ptr = nullptr;
            fixAlpha = 0;
            phredTable = getPhredTable();
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
            phredTable = getPhredTable();
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

        // Conditional base likelihood lookup table indexed by [is_error][genotype][base_class].
        // base_class: 0=ref (./,), 1=alt (matches altBase), 2=other
        static constexpr double COND_LK[2][3][3] = {
            // is_error = 0 (no sequencing error)
            {
                /* geno 0 (hom ref) */ {1.0,     0.0,     0.0},
                /* geno 1 (het)     */ {0.5,     0.5,     0.0},
                /* geno 2 (hom alt) */ {0.0,     1.0,     0.0},
            },
            // is_error = 1 (sequencing error)
            {
                /* geno 0 (hom ref) */ {0.0,     1.0/3.0, 2.0/3.0},
                /* geno 1 (het)     */ {1.0/6.0, 1.0/6.0, 2.0/3.0},
                /* geno 2 (hom alt) */ {1.0/3.0, 0.0,     2.0/3.0},
            },
        };

        // Classify a base as ref (0), alt (1), or other (2).
        static inline int classifyBase(char base, char altBase) {
            if (base == '.' || base == ',') return 0;
            if (toupper(base) == toupper(altBase)) return 1;
            return 2;
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

            // Precompute log-likelihood contributions for every possible combination
            // of (base_class, quality_score, g1, g2).
            //
            // The per-base log-likelihood under a genotype pair (g1, g2) is:
            //   log( [alpha*P(b|g1,err) + (1-alpha)*P(b|g2,err)] * P(err)
            //      + [alpha*P(b|g1,ok)  + (1-alpha)*P(b|g2,ok) ] * P(ok) )
            //
            // Within a single call alpha is constant, and the remaining inputs are
            // all discrete with small domains:
            //   - base_class (ref/alt/other):  3 values — determines P(b|g,err/ok)
            //   - quality score (Phred 0-93): 94 values — determines P(err), P(ok)
            //   - g1, g2 (genotype pair):      3 x 3 = 9 combinations
            //
            // Total: 3 x 94 x 9 = 2,538 entries.  Building this table costs 2,538
            // log() calls, but replaces millions (depth x markers x 9) in the inner
            // loop — typically a ~500x reduction in log() calls.
            const double oneMinusAlpha = 1.0 - alpha;
            double logLkTable[3][94][3][3];
            for (int bc = 0; bc < 3; ++bc) {
                const double lkErr[3]   = { COND_LK[1][0][bc], COND_LK[1][1][bc], COND_LK[1][2][bc] };
                const double lkNoErr[3] = { COND_LK[0][0][bc], COND_LK[0][1][bc], COND_LK[0][2][bc] };
                for (int q = 0; q < 94; ++q) {
                    const double pErr   = phredTable[q];
                    const double pNoErr = 1.0 - pErr;
                    for (int g1 = 0; g1 < 3; ++g1) {
                        for (int g2 = 0; g2 < 3; ++g2) {
                            double val = (alpha * lkErr[g1] + oneMinusAlpha * lkErr[g2]) * pErr
                                       + (alpha * lkNoErr[g1] + oneMinusAlpha * lkNoErr[g2]) * pNoErr;
                            logLkTable[bc][q][g1][g2] = log(val);
                        }
                    }
                }
            }

            double sumLLK(0);
#ifdef _OPENMP
            omp_set_num_threads(ptr->numThread);
#pragma omp parallel for reduction (+:sumLLK)
#endif
            for (size_t i = 0; i < ptr->NumMarker; ++i) {

                const ContaminationEstimator::ResolvedMarker& rm = ptr->resolvedMarkers[i];
                if (rm.baseInfoIndex < 0) continue;

                const std::vector<char>& tmpBase = ptr->viewer.baseInfo[rm.baseInfoIndex];
                const std::vector<char>& tmpQual = ptr->viewer.qualInfo[rm.baseInfoIndex];

                if (tmpBase.size() == 0) continue;

                if (not ptr->isSanityCheckDisabled and
                    (tmpBase.size() < (ptr->viewer.avgDepth - 3 * ptr->viewer.sdDepth) or
                     tmpBase.size() > (ptr->viewer.avgDepth + 3 * ptr->viewer.sdDepth)))
                    continue;

                if (ptr->isAFknown) {
                    ptr->AFs[i] = ptr->AF2s[i] = rm.knownAFValue;
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

                char altBase = rm.altBase;

                // Accumulate log-likelihoods across all observed bases at this marker.
                // For each base we classify it (ref/alt/other), extract its quality
                // score, and use those two values to index into the precomputed
                // logLkTable — one addition per genotype pair instead of a log() call.
                //
                // baseLKAccum[g1][g2] = sum_j log P(base_j | g1, g2, alpha)
                //                     = log product_j P(base_j | g1, g2, alpha)
                const int depth = tmpBase.size();
                double baseLKAccum[3][3] = {};

                for (int j = 0; j < depth; ++j) {
                    const int bc = classifyBase(tmpBase[j], altBase);
                    // Clamp quality to [0, 93] before indexing the lookup
                    // table.  BAM-derived qualities are guaranteed in-range,
                    // but --PileupFile inputs are ASCII and may carry values
                    // outside the expected Phred+33 window if the file is
                    // malformed or uses a different encoding; clamping avoids
                    // out-of-bounds access without rejecting the record.
                    int q = static_cast<unsigned char>(tmpQual[j]) - 33;
                    if (q < 0) q = 0;
                    else if (q > 93) q = 93;

                    for (int g1 = 0; g1 < 3; ++g1)
                        for (int g2 = 0; g2 < 3; ++g2)
                            baseLKAccum[g1][g2] += logLkTable[bc][q][g1][g2];
                }

                // Marginalize over genotype pairs, weighting each by its prior
                // probability under Hardy-Weinberg (GF for contaminating, GF2 for intended).
                for (int g1 = 0; g1 < 3; ++g1)
                    for (int g2 = 0; g2 < 3; ++g2)
                        markerLK += exp(baseLKAccum[g1][g2]) * GF[g1] * GF2[g2];
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
//                std::cerr << "ContaminatingSamplePC1:" << globalPC[0] << "\tContaminatingSamplePC2:" << globalPC[1]
//                          << "\tIntendedSamplePC1:" << globalPC2[0] << "\tIntendedSamplePC2:" << globalPC2[1]
//                          << "\tFREEMIX(Alpha):" << globalAlpha << "\tllk:" << llk1 << std::endl;
                notice("ContaminatingSamplePC1:%f\tContaminatingSamplePC2:%f\tIntendedSamplePC1:%f\tIntendedSamplePC2:%f\tFREEMIX(Alpha):%f\tllk:%f",
                       globalPC[0],globalPC[1],globalPC2[0],globalPC2[1],globalAlpha,llk1);
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

    // Per-marker data resolved once by BuildResolvedMarkers() so that
    // ComputeMixLLKs can use flat vector indexing instead of repeated
    // nested hash map lookups (posIndex, ChooseBed, knownAF).
    struct ResolvedMarker {
        int baseInfoIndex;   // index into viewer.baseInfo/qualInfo, or -1 if marker absent
        char altBase;        // alternate allele from ChooseBed
        double knownAFValue; // allele frequency from knownAF (0.0 if !isAFknown)
    };
    std::vector<ResolvedMarker> resolvedMarkers;

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

    int ReadBam(const char *bamFile, const char *faiFile, const char *bedFile,
                mplp_conf_t *mplpPtr);

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

    void BuildResolvedMarkers();

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
