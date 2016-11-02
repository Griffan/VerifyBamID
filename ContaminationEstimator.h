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
    bool isHeter;
    int numPC;
    int seed;
#define PCtype double

//#define PHRED(x)    pow(10.0,x/-10.0)
    static double Phred(double x) {
        return pow(10.0, x / -10.0);
    }

    class fullLLKFunc : public VectorFunc {
    public:
        double min_af;
        double max_af;
        double llk;
        ContaminationEstimator *ptr;
        vector<double> localPC;
        vector<double> localPC2;
        double localAlpha;
        const char* Base;
        fullLLKFunc()
        {
            fullLLKFunc::Base = "actg";
            min_af=0.0005;
            max_af=0.9995;
            llk=0;
            ptr=nullptr;
            localAlpha=0;
        }

        fullLLKFunc(int dim, ContaminationEstimator* contPtr):localPC(dim,0.),localPC2(dim,0.) {
            fullLLKFunc::Base = "actg";
            min_af = 0.0005;
            max_af = 0.9995;
            llk = 0;
            ptr = contPtr;
            localAlpha = 0;
        }

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

        void InitialGF(double AF, double *GF) const {
            if (AF < min_af) AF = min_af;
            if (AF > max_af) AF = max_af;
            GF[0] = (1 - AF) * (1 - AF);
            GF[1] = 2 * (AF) * (1 - AF);
            GF[2] = AF * AF;
        }

        inline double computeMixLLKs(const vector<double> & tPC1,const vector<double> & tPC2,const double alpha)
        {

            double sumLLK(0);
#ifdef _OPENMP
            omp_set_num_threads(16);
#pragma omp parallel for reduction (+:sumLLK)
#endif
            for (size_t i = 0; i < ptr->NumMarker; ++i) {
                double markerLK(0);
                double GF[3];
                double GF2[3];
                std::string chr = ptr->PosVec[i].first;
                int pos = ptr->PosVec[i].second;
                if (ptr->viewer.posIndex.find(chr) == ptr->viewer.posIndex.end()) {
                    continue;
                }
                else if (ptr->viewer.posIndex[chr].find(pos) == ptr->viewer.posIndex[chr].end()) {
                    continue;
                }

                for (int k = 0; k <tPC1.size(); ++k) {
                    ptr->AFs[i]+=ptr->UD[i][k] * tPC1[k];
                }
                ptr->AFs[i] += ptr->means[i];
                ptr->AFs[i] /= 2.0;

                for (int k = 0; k <tPC2.size(); ++k) {
                    ptr->AF2s[i]+=ptr->UD[i][k] * tPC2[k];
                }
                ptr->AF2s[i] += ptr->means[i];
                ptr->AF2s[i] /= 2.0;


                InitialGF(ptr->AFs[i], GF);
                InitialGF(ptr->AF2s[i], GF2);
                std::vector<char> tmpBase = ptr->viewer.GetBaseInfoAt(chr, pos);
                std::vector<char> tmpQual = ptr->viewer.GetQualInfoAt(chr, pos);
                if (tmpBase.size() == 0) continue;

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
//                            std::cerr <<i<<"th marker\t"<<tmpBase[j]<<"\t"<<tmpQual[j]<<"\t"<<altBase<<"\tlocalAlpha:"<<localAlpha<<"\tgeno1:"<<geno1<<"\tgeno2:"<<geno2
//                            <<"\tgetConditionalBaseLK1:"<<getConditionalBaseLK(tmpBase[j], geno1, altBase, 1)<<"\t"<< getConditionalBaseLK(tmpBase[j], geno2, altBase, 1)<<"\tPhred:"<<Phred(tmpQual[j] - 33)
//                            <<"\tgetConditionalBaseLK0:"<<getConditionalBaseLK(tmpBase[j], geno1, altBase, 0)<<"\t"<<getConditionalBaseLK(tmpBase[j], geno2, altBase, 0)<< std::endl;
                        }
//                        std::cerr<<"baseLK:"<<baseLK;
                        markerLK += exp(baseLK) * GF[geno1] * GF2[geno2];
//                        std::cerr<<"markerLK:"<<markerLK<<std::endl;
                    }
                if (markerLK > 0)
                    sumLLK += log(markerLK);
//                std::cerr << "sumLLK:" << sumLLK << std::endl;
            }
//            std::cerr << "sumLLK:" << sumLLK << std::endl;
            return sumLLK;
        }
        inline double computeMixLLKs(const vector<double> & tPC, const double alpha)
        {

            double sumLLK(0);
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

                for (int k = 0; k <tPC.size(); ++k) {
                    ptr->AFs[i]+=ptr->UD[i][k] * tPC[k];
                }
                ptr->AFs[i] += ptr->means[i];
                ptr->AFs[i] /= 2.0;


                InitialGF(ptr->AFs[i], GF);

                std::vector<char> tmpBase = ptr->viewer.GetBaseInfoAt(chr, pos);
                std::vector<char> tmpQual = ptr->viewer.GetQualInfoAt(chr, pos);
                if (tmpBase.size() == 0) continue;

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
//                            std::cerr <<i<<"th marker\t"<<tmpBase[j]<<"\t"<<tmpQual[j]<<"\t"<<altBase<<"\tlocalAlpha:"<<localAlpha<<"\tgeno1:"<<geno1<<"\tgeno2:"<<geno2
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

        inline double computeMixLLKs(const vector<double> & tPC1, const vector<double> & tPC2)//fix localAlpha
        {

            double sumLLK(0);

#ifdef _OPENMP
            omp_set_num_threads(16);
#pragma omp parallel for reduction (+:sumLLK)
#endif
            for (size_t i = 0; i < ptr->NumMarker; ++i) {
                double markerLK(0);
                double GF[3];
                double GF2[3];
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
                {
                    for (int k = 0; k <tPC1.size(); ++k) {
                        ptr->AFs[i]+=ptr->UD[i][k] * tPC1[k];
                    }
                    ptr->AFs[i] += ptr->means[i];
                    ptr->AFs[i] /= 2.0;
                    for (int k = 0; k <tPC2.size(); ++k) {
                        ptr->AF2s[i]+=ptr->UD[i][k] * tPC2[k];
                    }
                    ptr->AF2s[i] += ptr->means[i];
                    ptr->AF2s[i] /= 2.0;
                }

                InitialGF(ptr->AFs[i], GF);
                InitialGF(ptr->AF2s[i], GF2);

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
                            baseLK += log(((1. - localAlpha) * getConditionalBaseLK(tmpBase[j], geno1, altBase, 1) +
                                           localAlpha * getConditionalBaseLK(tmpBase[j], geno2, altBase, 1)) *
                                          Phred(tmpQual[j] - 33)
                                          + ((1. - localAlpha) * getConditionalBaseLK(tmpBase[j], geno1, altBase, 0) +
                                             localAlpha * getConditionalBaseLK(tmpBase[j], geno2, altBase, 0)) *
                                            (1 - Phred(tmpQual[j] - 33)));
//                            std::cerr <<i<<"th marker\t"<<tmpBase[j]<<"\t"<<tmpQual[j]<<"\t"<<altBase<<"\tlocalAlpha:"<<localAlpha<<"\tgeno1:"<<geno1<<"\tgeno2:"<<geno2
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
        inline double computeMixLLKs(const vector<double> & tPC)//fix localAlpha
        {

            double sumLLK(0);

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
                {
                    for (int k = 0; k <tPC.size(); ++k) {
                        ptr->AFs[i]+=ptr->UD[i][k] * tPC[k];
                    }
                    ptr->AFs[i] += ptr->means[i];
                    ptr->AFs[i] /= 2.0;
                }

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
                            baseLK += log(((1. - localAlpha) * getConditionalBaseLK(tmpBase[j], geno1, altBase, 1) +
                                           localAlpha * getConditionalBaseLK(tmpBase[j], geno2, altBase, 1)) *
                                          Phred(tmpQual[j] - 33)
                                          + ((1. - localAlpha) * getConditionalBaseLK(tmpBase[j], geno1, altBase, 0) +
                                             localAlpha * getConditionalBaseLK(tmpBase[j], geno2, altBase, 0)) *
                                            (1 - Phred(tmpQual[j] - 33)));
//                            std::cerr <<i<<"th marker\t"<<tmpBase[j]<<"\t"<<tmpQual[j]<<"\t"<<altBase<<"\tlocalAlpha:"<<localAlpha<<"\tgeno1:"<<geno1<<"\tgeno2:"<<geno2
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

        inline double computeMixLLKs_mix(const double alpha)//fix PCs for two populations mixture
        {

            double sumLLK(0);
#ifdef _OPENMP
            omp_set_num_threads(16);
#pragma omp parallel for reduction (+:sumLLK)
#endif
            for (size_t i = 0; i < ptr->NumMarker; ++i) {
                double markerLK(0);
                double GF[3];
                double GF2[3];
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
                {
                    for (int k = 0; k <localPC.size(); ++k) {
                        ptr->AFs[i]+=ptr->UD[i][k] * localPC[k];
                    }
                    ptr->AFs[i] += ptr->means[i];
                    ptr->AFs[i] /= 2.0;
                    for (int k = 0; k <localPC2.size(); ++k) {
                        ptr->AF2s[i]+=ptr->UD[i][k] * localPC2[k];
                    }
                    ptr->AF2s[i] += ptr->means[i];
                    ptr->AF2s[i] /= 2.0;
                }

                InitialGF(ptr->AFs[i], GF);
                InitialGF(ptr->AF2s[i], GF2);

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
//                            std::cerr <<i<<"th marker\t"<<tmpBase[j]<<"\t"<<tmpQual[j]<<"\t"<<altBase<<"\tlocalAlpha:"<<localAlpha<<"\tgeno1:"<<geno1<<"\tgeno2:"<<geno2
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
        inline double computeMixLLKs(const double alpha)//fix PCs
        {

            double sumLLK(0);
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
                {
                    for (int k = 0; k <localPC.size(); ++k) {
                        ptr->AFs[i]+=ptr->UD[i][k] * localPC[k];
                    }
                    ptr->AFs[i] += ptr->means[i];
                    ptr->AFs[i] /= 2.0;
                }


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
//                            std::cerr <<i<<"th marker\t"<<tmpBase[j]<<"\t"<<tmpQual[j]<<"\t"<<altBase<<"\tlocalAlpha:"<<localAlpha<<"\tgeno1:"<<geno1<<"\tgeno2:"<<geno2
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
//            srand(static_cast<unsigned>(time(NULL)));
            for (int k = 0; k <localPC.size(); ++k) {
                localPC[k]=static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            }
            for (int k = 0; k <localPC2.size(); ++k) {
                localPC2[k]=static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            }
            localAlpha = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);

            if(!ptr->isHeter)
            {
                if(ptr->isPCFixed) {
                    for (int k = 0; k <localPC.size(); ++k) {
                        localPC[k]=ptr->PC[0][k];
                    }
                    llk = (0 - computeMixLLKs(localAlpha));
                }
                else if(ptr->isAlphaFixed) {
                    localAlpha=ptr->alpha;
                    llk = (0 - computeMixLLKs(localPC));
                }
                else
                    llk = (0 - computeMixLLKs(localPC,localAlpha));
            }
            else//contamination source from different population
            {
                if(ptr->isPCFixed) {
                    for (int k = 0; k <localPC.size(); ++k) {
                        localPC[k]=ptr->PC[0][k];
                    }
                    for (int k = 0; k <localPC2.size(); ++k) {
                        localPC[k]=ptr->PC[1][k];
                    }
                    llk = (0 - computeMixLLKs_mix(localAlpha));
                }
                else if(ptr->isAlphaFixed) {
                    localAlpha=ptr->alpha;
                    llk = (0 - computeMixLLKs(localPC, localPC2));
                }
                else
                    llk = (0 - computeMixLLKs(localPC,localPC2,localAlpha));
            }
        }

        int initialize() {
            srand(ptr->seed);
//            srand(static_cast<unsigned>(time(NULL)));
            for (int k = 0; k <localPC.size(); ++k) {
                localPC[k]=static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            }
            for (int k = 0; k <localPC2.size(); ++k) {
                localPC2[k]=static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            }
            localAlpha = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);

            if(!ptr->isHeter)
            {
                if(ptr->isPCFixed) {
                    for (int k = 0; k <localPC.size(); ++k) {
                        localPC[k]=ptr->PC[0][k];
                    }
                    llk = (0 - computeMixLLKs(localAlpha));
                }
                else if(ptr->isAlphaFixed) {
                    localAlpha=ptr->alpha;
                    llk = (0 - computeMixLLKs(localPC));
                }
                else
                    llk = (0 - computeMixLLKs(localPC,localAlpha));
            }
            else//contamination source from different population
            {
                if(ptr->isPCFixed) {
                    for (int k = 0; k <localPC.size(); ++k) {
                        localPC[k]=ptr->PC[0][k];
                    }
                    for (int k = 0; k <localPC2.size(); ++k) {
                        localPC[k]=ptr->PC[1][k];
                    }
                    llk = (0 - computeMixLLKs_mix(localAlpha));
                }
                else if(ptr->isAlphaFixed) {
                    localAlpha=ptr->alpha;
                    llk = (0 - computeMixLLKs(localPC, localPC2));
                }
                else
                    llk = (0 - computeMixLLKs(localPC,localPC2,localAlpha));
            }
            return 0;
        }

        virtual double Evaluate(Vector &v) {
            double smLLK = 0;
            if(!ptr->isHeter)
            {
                if(ptr->isPCFixed) {
                    double tmpAlpha = invLogit(v[0]);
                    smLLK = 0 - computeMixLLKs(tmpAlpha);
                    if (smLLK < llk) {
                        llk = smLLK;
                        localAlpha = tmpAlpha;
                    }
                }
                else if(ptr->isAlphaFixed) {
                    vector<double> tmpPC(ptr->numPC,0.);
                    for (int i = 0; i <ptr->numPC ; ++i) {
                        tmpPC[i]=v[i];
                    }
                    smLLK = 0 - computeMixLLKs(tmpPC);
                    if (smLLK < llk) {
                        llk = smLLK;
                        localPC = tmpPC;
                    }
                }
                else {
                    vector<double> tmpPC(ptr->numPC,0.);
                    for (int i = 0; i <ptr->numPC ; ++i) {
                        tmpPC[i]=v[i];
                    }
                    double tmpAlpha=invLogit(v[ptr->numPC]);
                    smLLK = (0 - computeMixLLKs(tmpPC, tmpAlpha));
                    if (smLLK < llk) {
                        llk = smLLK;
                        localPC = tmpPC;
                        localAlpha = tmpAlpha;
                    }
                }
            }
            else//contamination source from different population
            {
                if(ptr->isPCFixed) {
                    double tmpAlpha = invLogit(v[0]);
                    smLLK = (0 - computeMixLLKs_mix(localAlpha));
                    if (smLLK < llk) {
                        llk = smLLK;
                        localAlpha = tmpAlpha;
                    }
                }
                else if(ptr->isAlphaFixed) {
                    vector<double> tmpPC(ptr->numPC,0.);
                    vector<double> tmpPC2(ptr->numPC,0.);

                    for (int k = 0; k <v.Length(); ++k) {
                        if(k < ptr->numPC)
                            tmpPC[k]=v[k];
                        else if(k< ptr->numPC*2)
                            tmpPC2[k-(ptr->numPC)]=v[k];
                        else {
                            error("Simplex Vector dimension error!");
                            exit(EXIT_FAILURE);
                        }
                    }
                    smLLK = (0 - computeMixLLKs(tmpPC, tmpPC2));
                    if (smLLK < llk) {
                        llk = smLLK;
                        localPC = tmpPC;
                        localPC2 = tmpPC2;
                    }
                }
                else {
                    vector<double> tmpPC(ptr->numPC,0.);
                    vector<double> tmpPC2(ptr->numPC,0.);
                    double tmpAlpha(0.);
                    for (int k = 0; k <v.Length(); ++k) {
                        if(k < ptr->numPC)
                            tmpPC[k]=v[k];
                        else if(k < ptr->numPC*2)
                            tmpPC2[k-(ptr->numPC)]=v[k];
                        else if(k == ptr->numPC*2)
                            tmpAlpha = invLogit(v[k]);
                        else{
                            error("Simplex Vector dimension error!");
                            exit(EXIT_FAILURE);
                        }
                    }
                    smLLK = (0 - computeMixLLKs(tmpPC, tmpPC2, tmpAlpha));
                    if (smLLK < llk) {
                        llk = smLLK;
                        localPC = tmpPC;
                        localPC2 = tmpPC2;
                        localAlpha = tmpAlpha;
                    }
                }
            }
            std::cerr << "tmpPC1:" << localPC[0] << "\ttmpPC2:" << localPC[1]
                      << "\ttmpPC3:" << localPC2[0] << "\ttmpPC4:" << localPC2[1]
                      << "\talpha:" << localAlpha << "\tllk:" << llk << std::endl;
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
    std::vector<double> AF2s;

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
    void OptimizeHomFixedPC(AmoebaMinimizer &myMinimizer);

    void OptimizeHomFixedAlpha(AmoebaMinimizer &myMinimizer);

    void OptimizeHom(AmoebaMinimizer &myMinimizer);

    void OptimizeHeterFixedPC(AmoebaMinimizer &myMinimizer);

    void OptimizeHeterFixedAlpha(AmoebaMinimizer &myMinimizer);

    void OptimizeHeter(AmoebaMinimizer &myMinimizer);
};

#endif /* CONTAMINATIONESTIMATOR_H_ */
