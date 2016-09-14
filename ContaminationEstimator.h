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

/* Contact: Fan Zhang <fanzhang@umich.edu> */
#ifndef CONTAMINATIONESTIMATOR_H_
#define CONTAMINATIONESTIMATOR_H_

#include <string>

#include <unordered_map>
//#include <tkDecls.h>
#include "MathVector.h"
#include "MathGenMin.h"
#include "SimplePileupViewer.h"
#ifdef _OPENMP
#include "omp.h"
#endif

class ContaminationEstimator {
public:
    bool isPCFixed;
    bool isAlphaFixed;
    bool isAFknown;
    int seed;
#define PCtype double

//#define PHRED(x)    pow(10.0,x/-10.0)
    static double Phred(double x) {
        return pow(10.0, x / -10.0);
    }

    class fullLLKFunc : public VectorFunc {
    public:

        double llk;
        ContaminationEstimator *ptr;
        double PC1, PC2, PC3, PC4,alpha;

        fullLLKFunc() { };

        ~fullLLKFunc() { };

        inline static double invLogit(double &x) {
            double e = exp(x);
            return e / (1. + e);
        };
/*
        inline double computeMixLLKs(double tPC1, double tPC2) {
            double min_af(0.5 / ptr->NumIndividual), max_af((ptr->NumIndividual - 0.5) / ptr->NumIndividual);
            double sumLLK(0), GF0(0), GF1(0), GF2(0);
            std::string chr;
            uint32_t pos;
            size_t glIndex = 0;
            for (size_t i = 0; i != ptr->NumMarker; ++i) {
                //std::cerr << "Number " << i << "th marker out of " << ptr->NumMarker << " markers and " << ptr->NumIndividual << " individuals"<<std::endl;
                //std::cerr << "AF:" << ptr->AFs[i] << "\tUD:" << ptr->UD[i][0] << "\t" << ptr->UD[i][1] << "\tmeans:" << ptr->means[i] << std::endl;
                chr = ptr->PosVec[i].first;
                pos = ptr->PosVec[i].second;
                if (ptr->MarkerIndex[chr].find(pos) != ptr->MarkerIndex[chr].end())
                    glIndex = ptr->MarkerIndex[chr][pos];
                else
                    glIndex = ptr->GL.size() - 1;
                ptr->AFs[i] = ((ptr->UD[i][0] * tPC1 + ptr->UD[i][1] * tPC2) + ptr->means[i]) / 2.0;
                if (ptr->AFs[i] < min_af) ptr->AFs[i] = min_af;
                if (ptr->AFs[i] > max_af) ptr->AFs[i] = max_af;
                GF0 = (1 - ptr->AFs[i]) * (1 - ptr->AFs[i]);
                GF1 = 2 * (ptr->AFs[i]) * (1 - ptr->AFs[i]);
                GF2 = (ptr->AFs[i]) * (ptr->AFs[i]);
                sumLLK += log(Phred(ptr->GL[glIndex][0]) * GF0 + Phred(ptr->GL[glIndex][1]) * GF1 +
                              Phred(ptr->GL[glIndex][2]) * GF2);
                //std::cerr << "GL:" << ptr->GL[glIndex][0] << "\t" << ptr->GL[glIndex][1] << "\t" << ptr->GL[glIndex][2] << "\t" << chr << "\t" << pos << std::endl;
                //std::cerr << "AF:" << ptr->AFs[i] << "\tUD:" << ptr->UD[i][0] << "\t" << ptr->UD[i][1] << "\tmeans:" << ptr->means[i] << std::endl;
            }
            //std::cerr << "sumLLK:" << sumLLK << std::endl;
            return sumLLK;
        }
*/
        static const std::string base = "actg";

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
            return base[maxIndex];
        }

        inline double getConditionalBaseLK(char base, int genotype, char altBase, bool is_error) {
            if (!is_error) {
                if (genotype == 0) {
                    if (base == '.' || base == ',') {
                        return 1;
                    }
                    else
                        return 0;
                }
                else if (genotype == 1) {
                    if (base == '.' || base == ',') {
                        return 0.5;
                    }
                    else if (toupper(base) == toupper(altBase)) {
                        return 0.5;
                    }
                    else
                        return 0;
                }
                else if (genotype == 2) {
                    if (toupper(base) == toupper(altBase)) {
                        return 1;
                    }
                    else
                        return 0;
                }
                else {
                    std::cerr << "genotype error!" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                if (genotype == 0) {
                    if (base == '.' || base == ',') {
                        return 0;
                    }
                    else if (toupper(base) == toupper(altBase)) {
                        return 1. / 3.;
                    }
                    else
                        return 2. / 3.;
                }
                else if (genotype == 1) {
                    if (base == '.' || base == ',') {
                        return 1. / 6.;
                    }
                    else if (toupper(base) == toupper(altBase)) {
                        return 1. / 6.;
                    }
                    else
                        return 2. / 3.;
                }
                else if (genotype == 2) {
                    if (base == '.' || base == ',') {
                        return 1. / 3.;
                    }
                    if (toupper(base) == toupper(altBase)) {
                        return 0;
                    }
                    else
                        return 2. / 3.;
                }
                else {
                    std::cerr << "genotype error!" << std::endl;
                    exit(EXIT_FAILURE);
                }

            }


        }

        inline double computeMixLLKs(double tPC1, double tPC2, double tPC3, double tPC4, double alpha) {
            double min_af(0.0005), max_af(0.9995);
            double sumLLK(0);//, GF0(0), GF1(0), GF2(0);
//            size_t glIndex = 0;
#ifdef _OPENMP
            omp_set_num_threads(16);
#pragma omp parallel for reduction (+:sumLLK)
#endif
            for (size_t i = 0; i < ptr->NumMarker; ++i) {
                double markerLK(0);
                double GF[3];
                std::string chr = ptr->PosVec[i].first;
                int pos = ptr->PosVec[i].second;
                if (ptr->viewer.posIndex.find(chr) == ptr->viewer.posIndex.end()) {
                    continue;
                }
                else if (ptr->viewer.posIndex[chr].find(pos) == ptr->viewer.posIndex[chr].end()) {
                    continue;
                }

                ptr->AFs[i] = ((ptr->UD[i][0] * tPC1 + ptr->UD[i][1] * tPC2+ ptr->UD[i][2] * tPC3 + ptr->UD[i][3] * tPC4) + ptr->means[i]) / 2.0;
//                std::cerr << "AF:" << ptr->AFs[i] << "\tUD:" << ptr->UD[i][0] << "\t" << ptr->UD[i][1] << "PC:"<<tPC1<<"\t"<<tPC2<<"\t"<<tPC3<<"\t"<<tPC4<<std::endl;

                if (ptr->AFs[i] < min_af) ptr->AFs[i] = min_af;
                if (ptr->AFs[i] > max_af) ptr->AFs[i] = max_af;
                GF[0] = (1 - ptr->AFs[i]) * (1 - ptr->AFs[i]);//genotype frequency
                GF[1] = 2 * (ptr->AFs[i]) * (1 - ptr->AFs[i]);
                GF[2] = (ptr->AFs[i]) * (ptr->AFs[i]);
                std::vector<char> tmpBase = ptr->viewer.GetBaseInfoAt(chr, pos);
                std::vector<char> tmpQual = ptr->viewer.GetQualInfoAt(chr, pos);
                if (tmpBase.size() == 0) continue;
//                std::cerr<<"chr:"<<chr<<"\tpos:"<<pos<<"\t";
//                std::cerr<<"tmpBase:";
//                for (int k = 0; k <tmpBase.size() ; ++k) {
//                    std::cerr<<tmpBase[k];
//                }
//                std::cerr<<"\t";
                //char altBase=findAlt(tmpBase);
                char altBase = ptr->ChooseBed[chr][pos].second;

                for (int geno1 = 0; geno1 < 3; ++geno1)
                    for (int geno2 = 0; geno2 < 3; ++geno2) {
                        double baseLK(0);
                        for (int j = 0; j < tmpBase.size(); ++j) {
                            baseLK += log(((1. - alpha) * getConditionalBaseLK(tmpBase[j], geno1, altBase, 1) +
                                           alpha * getConditionalBaseLK(tmpBase[j], geno2, altBase, 1)) *
                                          Phred(tmpQual[j] - 33)
                                          + ((1. - alpha) * getConditionalBaseLK(tmpBase[j], geno1, altBase, 0) +
                                             alpha * getConditionalBaseLK(tmpBase[j], geno2, altBase, 0)) *
                                            (1 - Phred(tmpQual[j] - 33)));
//                            std::cerr <<i<<"th marker\t"<<tmpBase[j]<<"\t"<<tmpQual[j]<<"\t"<<altBase<<"\talpha:"<<alpha<<"\tgeno1:"<<geno1<<"\tgeno2:"<<geno2
//                            <<"\tgetConditionalBaseLK1:"<<getConditionalBaseLK(tmpBase[j], geno1, altBase, 1)<<"\t"<< getConditionalBaseLK(tmpBase[j], geno2, altBase, 1)<<"\tPhred:"<<Phred(tmpQual[j] - 33)
//                            <<"\tgetConditionalBaseLK0:"<<getConditionalBaseLK(tmpBase[j], geno1, altBase, 0)<<"\t"<<getConditionalBaseLK(tmpBase[j], geno2, altBase, 0)<< std::endl;
                        }
//                        std::cerr<<"baseLK:"<<baseLK;
                        markerLK += exp(baseLK) * GF[geno1] * GF[geno2];
//                        std::cerr<<"markerLK:"<<markerLK<<std::endl;
                    }
                if (markerLK > 0)
                    sumLLK += log(markerLK);
//                std::cerr << "sumLLK:" << sumLLK << std::endl;
            }
//            std::cerr << "sumLLK:" << sumLLK << std::endl;
            return sumLLK;
        }

        inline double computeMixLLKs(double tPC1, double tPC2) {
            double min_af(0.0005), max_af(0.9995);
            alpha=0.01;
            double sumLLK(0);//, GF0(0), GF1(0), GF2(0);
//            size_t glIndex = 0;
#ifdef _OPENMP
            omp_set_num_threads(16);
#pragma omp parallel for reduction (+:sumLLK)
#endif
            for (size_t i = 0; i < ptr->NumMarker; ++i) {
                double markerLK(0);
                double GF[3];
                std::string chr = ptr->PosVec[i].first;
                int pos = ptr->PosVec[i].second;
                if (ptr->viewer.posIndex.find(chr) == ptr->viewer.posIndex.end()) {
                    continue;
                }
                else if (ptr->viewer.posIndex[chr].find(pos) == ptr->viewer.posIndex[chr].end()) {
                    continue;
                }
                if(ptr->isAFknown)
                {
                    ptr->AFs[i]=ptr->knownAF[chr][pos];
                }
                else
                    ptr->AFs[i] = ((ptr->UD[i][0] * tPC1 + ptr->UD[i][1] * tPC2) + ptr->means[i]) / 2.0;
//                std::cerr << "AF:" << ptr->AFs[i] << "\tUD:" << ptr->UD[i][0] << "\t" << ptr->UD[i][1] << "PC:"<<tPC1<<"\t"<<tPC2<<"\tmeans:" << ptr->means[i] << std::endl;

                if (ptr->AFs[i] < min_af) ptr->AFs[i] = min_af;
                if (ptr->AFs[i] > max_af) ptr->AFs[i] = max_af;
                GF[0] = (1 - ptr->AFs[i]) * (1 - ptr->AFs[i]);//genotype frequency
                GF[1] = 2 * (ptr->AFs[i]) * (1 - ptr->AFs[i]);
                GF[2] = (ptr->AFs[i]) * (ptr->AFs[i]);
                std::vector<char> tmpBase = ptr->viewer.GetBaseInfoAt(chr, pos);
                std::vector<char> tmpQual = ptr->viewer.GetQualInfoAt(chr, pos);
                if (tmpBase.size() == 0) continue;
//                std::cerr<<"chr:"<<chr<<"\tpos:"<<pos<<"\t";
//                std::cerr<<"tmpBase:";
//                for (int k = 0; k <tmpBase.size() ; ++k) {
//                    std::cerr<<tmpBase[k];
//                }
//                std::cerr<<"\t";
                //char altBase=findAlt(tmpBase);
                char altBase = ptr->ChooseBed[chr][pos].second;

                for (int geno1 = 0; geno1 < 3; ++geno1)
                    for (int geno2 = 0; geno2 < 3; ++geno2) {
                        double baseLK(0);
                        for (int j = 0; j < tmpBase.size(); ++j) {
                            baseLK += log(((1. - alpha) * getConditionalBaseLK(tmpBase[j], geno1, altBase, 1) +
                                           alpha * getConditionalBaseLK(tmpBase[j], geno2, altBase, 1)) *
                                          Phred(tmpQual[j] - 33)
                                          + ((1. - alpha) * getConditionalBaseLK(tmpBase[j], geno1, altBase, 0) +
                                             alpha * getConditionalBaseLK(tmpBase[j], geno2, altBase, 0)) *
                                            (1 - Phred(tmpQual[j] - 33)));
//                            std::cerr <<i<<"th marker\t"<<tmpBase[j]<<"\t"<<tmpQual[j]<<"\t"<<altBase<<"\talpha:"<<alpha<<"\tgeno1:"<<geno1<<"\tgeno2:"<<geno2
//                            <<"\tgetConditionalBaseLK1:"<<getConditionalBaseLK(tmpBase[j], geno1, altBase, 1)<<"\t"<< getConditionalBaseLK(tmpBase[j], geno2, altBase, 1)<<"\tPhred:"<<Phred(tmpQual[j] - 33)
//                            <<"\tgetConditionalBaseLK0:"<<getConditionalBaseLK(tmpBase[j], geno1, altBase, 0)<<"\t"<<getConditionalBaseLK(tmpBase[j], geno2, altBase, 0)<< std::endl;
                        }
//                        std::cerr<<"baseLK:"<<baseLK;
                        markerLK += exp(baseLK) * GF[geno1] * GF[geno2];
//                        std::cerr<<"markerLK:"<<markerLK<<std::endl;
                    }
                if (markerLK > 0)
                    sumLLK += log(markerLK);
//                std::cerr << "sumLLK:" << sumLLK << std::endl;
            }
//            std::cerr << "sumLLK:" << sumLLK << std::endl;
            return sumLLK;
        }

        inline double computeMixLLKs(double alpha) {
            double min_af(0.0005), max_af(0.9995);
            double sumLLK(0);//, GF0(0), GF1(0), GF2(0);
#ifdef _OPENMP
            omp_set_num_threads(16);
#pragma omp parallel for reduction (+:sumLLK)
#endif
            for (size_t i = 0; i < ptr->NumMarker; ++i) {
                double markerLK(0);
                double GF[3];
                std::string chr = ptr->PosVec[i].first;
                int pos = ptr->PosVec[i].second;
                if (ptr->viewer.posIndex.find(chr) == ptr->viewer.posIndex.end()) {
                    continue;
                }
                else if (ptr->viewer.posIndex[chr].find(pos) == ptr->viewer.posIndex[chr].end()) {
                    continue;
                }

                if(ptr->isAFknown)
                {
                    ptr->AFs[i]=ptr->knownAF[chr][pos];
                }
                else
                    ptr->AFs[i] = (ptr->means[i]) / 2.0;
//                std::cerr << "AF:" << ptr->AFs[i] <<"\tmeans:" << ptr->means[i]/2.0 << std::endl;

                if (ptr->AFs[i] < min_af) ptr->AFs[i] = min_af;
                if (ptr->AFs[i] > max_af) ptr->AFs[i] = max_af;
                GF[0] = (1 - ptr->AFs[i]) * (1 - ptr->AFs[i]);//genotype frequency
                GF[1] = 2 * (ptr->AFs[i]) * (1 - ptr->AFs[i]);
                GF[2] = (ptr->AFs[i]) * (ptr->AFs[i]);
                std::vector<char> tmpBase = ptr->viewer.GetBaseInfoAt(chr, pos);
                std::vector<char> tmpQual = ptr->viewer.GetQualInfoAt(chr, pos);
                if (tmpBase.size() == 0) continue;
//                std::cerr<<"chr:"<<chr<<"\tpos:"<<pos<<"\t";
//                std::cerr<<"tmpBase:";
//                for (int k = 0; k <tmpBase.size() ; ++k) {
//                    std::cerr<<tmpBase[k];
//                }
//                std::cerr<<"\t";
                //char altBase=findAlt(tmpBase);
                char altBase = ptr->ChooseBed[chr][pos].second;

                for (int geno1 = 0; geno1 < 3; ++geno1)
                    for (int geno2 = 0; geno2 < 3; ++geno2) {
                        double baseLK(0);
                        for (int j = 0; j < tmpBase.size(); ++j) {
                            baseLK += log(((1. - alpha) * getConditionalBaseLK(tmpBase[j], geno1, altBase, 1) +
                                           alpha * getConditionalBaseLK(tmpBase[j], geno2, altBase, 1)) *
                                          Phred(tmpQual[j] - 33)
                                          + ((1. - alpha) * getConditionalBaseLK(tmpBase[j], geno1, altBase, 0) +
                                             alpha * getConditionalBaseLK(tmpBase[j], geno2, altBase, 0)) *
                                            (1 - Phred(tmpQual[j] - 33)));
//                            std::cerr <<i<<"th marker\t"<<tmpBase[j]<<"\t"<<tmpQual[j]<<"\t"<<altBase<<"\talpha:"<<alpha<<"\tgeno1:"<<geno1<<"\tgeno2:"<<geno2
//                            <<"\tgetConditionalBaseLK1:"<<getConditionalBaseLK(tmpBase[j], geno1, altBase, 1)<<"\t"<< getConditionalBaseLK(tmpBase[j], geno2, altBase, 1)<<"\tPhred:"<<Phred(tmpQual[j] - 33)
//                            <<"\tgetConditionalBaseLK0:"<<getConditionalBaseLK(tmpBase[j], geno1, altBase, 0)<<"\t"<<getConditionalBaseLK(tmpBase[j], geno2, altBase, 0)<< std::endl;
                        }
//                        std::cerr<<"baseLK:"<<baseLK;
                        markerLK += exp(baseLK) * GF[geno1] * GF[geno2];
//                        std::cerr<<"markerLK:"<<markerLK<<std::endl;
                    }
                if (markerLK > 0)
                    sumLLK += log(markerLK);
//                std::cerr << "sumLLK:" << sumLLK << std::endl;
            }
//            std::cerr << "sumLLK:" << sumLLK << std::endl;
            return sumLLK;
        }

        fullLLKFunc(ContaminationEstimator *inPtr) {
            ptr = inPtr;
            srand(ptr->seed);
            double r1 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            double r2 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            double r3 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            double r4 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            double r5 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            if(ptr->isPCFixed)
                llk = (0 - computeMixLLKs(r3));
            else if(ptr->isAlphaFixed)
                llk = (0 - computeMixLLKs(r1,r2));
            else
                llk = (0 - computeMixLLKs(r1,r2,r3,r4,r5));
        }

        int initialize(ContaminationEstimator *inPtr) {
            ptr = inPtr;
            srand(ptr->seed);
            double r1 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            double r2 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            double r3 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            double r4 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            double r5 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            if(ptr->isPCFixed)
                llk = (0 - computeMixLLKs(r3));
            else if(ptr->isAlphaFixed)
                llk = (0 - computeMixLLKs(r1,r2));
            else
                llk = (0 - computeMixLLKs(r1,r2,r3,r4,r5));
            return 0;
        }

        virtual double Evaluate(Vector &v) {
            double smLLK = 0;
            switch (v.Length()) {
                case 5: {
                    double tmpPC1 = v[0];
                    double tmpPC2 = v[1];
                    double tmpPC3 = v[2];
                    double tmpPC4 = v[3];
                    double tmpAlpha = invLogit(v[4]);
                    smLLK = 0 - computeMixLLKs(tmpPC1, tmpPC2, tmpPC3, tmpPC4,tmpAlpha);
                    if (smLLK < llk) {
                        llk = smLLK;
                        PC1 = tmpPC1;
                        PC2 = tmpPC2;
                        PC3 = tmpPC3;
                        PC4 = tmpPC4;
                        alpha = tmpAlpha;
                    }
                    break;
                }

                case 2: {
                    double tmpPC1 = v[0];
                    double tmpPC2 = v[1];
                    smLLK = 0 - computeMixLLKs(tmpPC1, tmpPC2);
                    if (smLLK < llk) {
                        llk = smLLK;
                        PC1 = tmpPC1;
                        PC2 = tmpPC2;
                    }
                    break;
                }

                case 1: {
                    double tmpAlpha = invLogit(v[0]);;
                    smLLK = 0 - computeMixLLKs(tmpAlpha);
                    if (smLLK < llk) {
                        llk = smLLK;
                        alpha = tmpAlpha;
                    }
                    break;
                }

                default:
                    error("Simplex Vector dimension error!");
                    exit(EXIT_FAILURE);
            }
            std::cerr << "tmpPC1:" << PC1 << "\ttmpPC2:" << PC2 << "\talpha:" << alpha << "\tllk:" << llk << std::endl;
            return smLLK;
        }
    };

    SimplePileupViewer viewer;
    double alpha;
    uint32_t NumMarker;
    //uint32_t NumIndividual;
    fullLLKFunc fn;

    //std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t> > MarkerIndex;
    std::unordered_map<std::string, std::unordered_map<uint32_t, double> > knownAF;
    //std::unordered_map<std::string, uint32_t> IndividualIndex;

    std::vector<std::vector<PCtype> > UD;
    std::vector<std::vector<PCtype> > PC;
    //std::vector<std::vector<double> > GL;
    std::vector<PCtype> means;
    std::vector<double> AFs;

    std::unordered_map<std::string, std::unordered_map<int, std::pair<char, char> > > ChooseBed;
    std::vector<region_t> BedVec;
    std::vector<std::pair<std::string, int> > PosVec;

    ContaminationEstimator();

    ContaminationEstimator(const char *bamFile, const char *faiFile, const char *bedFile);

    /*Initialize from existed UD*/
    /*This assumes the markers are the same as the selected vcf*/
    /*ContaminationEstimator(const std::string &UDpath, const std::string &PCpath, const std::string &Mean,
                           const std::string &pileup, const std::string &GLpath, const std::string &Bed);
    */
    int ReadMatrixUD(const std::string &path);
    /*
    int ReadMatrixPC(const std::string &path);
    */
    /*Intersect marker sites*/
    /*
    int ReadMatrixGL(const std::string &path);
    */
    int ReadChooseBed(const std::string &path);

    int ReadMean(const std::string &path);

    int ReadAF(const std::string & path);
    /*
    int CheckMarkerSetConsistency();

    int FormatMarkerIntersection();
    */
    /*Optimize*/
    int OptimizeLLK();
    /*
    int RunMapping();

    int writeVcfFile(const std::string &path);

    int ReadPileup(const std::string &path);
    */
    ~ContaminationEstimator();
    /*
    int RunFromVCF(const std::string VcfSiteAFFile, const std::string CurrentMPU, const std::string ReadGroup,
                   const std::string Prefix);

    int RunFromSVDMatrix(const std::string UDpath, const std::string PCpath, const std::string Mean,
                         const std::string &MPUpath, const std::string &Bed, const std::string &Prefix,
                         const std::string &ReadGroup);
    */
    int ReadSVDMatrix(const std::string UDpath, const std::string Mean, const std::string &Bed);
    /*
    int FromBamtoPileup();
     */
};

#endif /* CONTAMINATIONESTIMATOR_H_ */
