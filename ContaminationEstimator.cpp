/*
 * ContaminationEstimator.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: fanzhang
 */

#include "ContaminationEstimator.h"
#include <fstream>
#include <sstream>
#include <cmath>

#ifndef MPU_PATH
#define MPU_PATH "mpuTools"
#endif

ContaminationEstimator::ContaminationEstimator() {

}

ContaminationEstimator::~ContaminationEstimator() {
}


int ContaminationEstimator::OptimizeLLK(const std::string &OutputPrefix) {
    AmoebaMinimizer myMinimizer;
    std::ofstream fout(OutputPrefix+".out");
    fn.initialize();
    std::vector<double> candidateAlphaSet{0.01,0.02,0.05,0.1,0.2};
    if (!isHeter) {
        if (isPCFixed) {
            std::cout << "Estimation from OptimizeHomoFixedPC:" << std::endl;
            fout<< "Estimation from OptimizeHomoFixedPC:" << std::endl;
            for (auto alphaTrial: candidateAlphaSet) {
                alpha = alphaTrial;//all the initializations in this function are for AmoebaMinimizer simplex
                if(OptimizeHomoFixedPC(myMinimizer)) break;
            }
        } else if (isAlphaFixed) {
            std::cout << "Estimation from OptimizeHomoFixedAlpha:" << std::endl;
            fout<< "Estimation from OptimizeHomoFixedAlpha:" << std::endl;
            for (int k = 0; k < PC[0].size(); ++k) {
                PC[0][k] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            }
            OptimizeHomoFixedAlpha(myMinimizer);
        } else {
            for (auto alphaTrial: candidateAlphaSet) {
                for (int k = 0; k < PC[0].size(); ++k) {
                    PC[0][k] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
                }
                alpha = alphaTrial;
                std::cout << "Estimation from OptimizeHomo:" << std::endl;
                fout << "Estimation from OptimizeHomo:" << std::endl;
                if(OptimizeHomo(myMinimizer)) break;
            }
        }
    } else//contamination source from different population
    {
        if (isPCFixed) {
            std::cout << "Estimation from OptimizeHeterFixedPC:" << std::endl;
            fout<< "Estimation from OptimizeHeterFixedPC:" << std::endl;
            for (auto alphaTrial: candidateAlphaSet) {
                alpha = alphaTrial;
                if(OptimizeHeterFixedPC(myMinimizer)) break;
            }
        } else if (isAlphaFixed) {
            for (int k = 0; k < numPC; ++k) {
                PC[0][k] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            }
            for (int k = 0; k < numPC; ++k) {
                PC[1][k] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            }
            std::cout << "Estimation from OptimizeHeterFixedAlpha:" << std::endl;
            fout<< "Estimation from OptimizeHeterFixedAlpha:" << std::endl;
            isHeter = false;
            OptimizeHomoFixedAlpha(myMinimizer);
            PC[1]=PC[0];
            fn.globalPC2 = fn.globalPC;
            isHeter = true;
            OptimizeHeterFixedAlpha(myMinimizer);
        } else {
            for (auto alphaTrial: candidateAlphaSet) {
                for (int k = 0; k < numPC; ++k) {
                    PC[0][k] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
                }
                alpha = alphaTrial;
                std::cout << "Estimation from OptimizeHeter:" << std::endl;
                fout << "Estimation from OptimizeHeter:" << std::endl;
                isHeter = false;
                OptimizeHomo(myMinimizer);
                PC[1] = PC[0];
                fn.globalPC2 = fn.globalPC;
                isHeter = true;
                if (OptimizeHeter(myMinimizer)) break;
            }
        }
	    if(fn.globalAlpha >= 0.5)
	    {
		    std::swap(fn.globalPC[0],fn.globalPC2[0]);
		    std::swap(fn.globalPC[1],fn.globalPC2[1]);
	    }
    }

    std::cout << "Contaminating Sample ";
    fout<< "Contaminating Sample ";
    for(int i =0; i < numPC; ++i)
    {
        std::cout <<"PC"<<i+1<<":" << fn.globalPC[i]<<"\t";
        fout<< "PC"<<i+1<<":" << fn.globalPC[i]<<"\t";
    }
    std::cout <<std::endl;
    fout<<std::endl;

    std::cout << "Intended Sample ";
    fout<< "Intended Sample ";
    for(int i =0; i < numPC; ++i)
    {
        std::cout <<"PC"<<i+1<<":" << fn.globalPC2[i]<<"\t";
        fout<<"PC"<<i+1<<":" << fn.globalPC2[i]<<"\t";
    }
    std::cout<< std::endl;
    fout<<std::endl;

    std::cout << "Alpha:" << (fn.globalAlpha < 0.5 ? fn.globalAlpha : (1 - fn.globalAlpha)) << std::endl;
    fout<< "Alpha:" << (fn.globalAlpha < 0.5 ? fn.globalAlpha : (1 - fn.globalAlpha)) << std::endl;
    fout.close();
    return 0;
}

bool ContaminationEstimator::OptimizeHeter(AmoebaMinimizer &myMinimizer) {
    Vector startingPoint("TestPoint", numPC * 2 + 1);
    for (int i = 0; i < numPC * 2; ++i) {
        if (i < numPC)
            startingPoint[i] = PC[0][i];
        else
            startingPoint[i] = PC[1][i - numPC];
    }
    startingPoint[numPC * 2] = fullLLKFunc::Logit(alpha);

    if(verbose) {
        std::cerr << "Start point:";
        for (int i = 0; i < numPC * 2; ++i) {
            std::cerr << startingPoint[i] << "\t";
        }
        std::cerr << "and alpha:\t" << alpha << std::endl;
    }

    startingPoint.label = "startPoint";
    myMinimizer.func = &fn;
    myMinimizer.Reset(numPC * 2 + 1);
    myMinimizer.point = startingPoint;
    double ret = myMinimizer.Minimize(epsilon);
    alpha = fullLLKFunc::invLogit(myMinimizer.point[numPC * 2]);
    for (int i = 0; i < numPC; ++i) {
        PC[0][i] = myMinimizer.point[i];
    }
    for (int i = numPC; i < numPC * 2; ++i) {
        PC[1][i - numPC] = myMinimizer.point[i];
    }
    if(ret == std::numeric_limits<double>::max()) return false;
    else return true;
}

bool ContaminationEstimator::OptimizeHeterFixedAlpha(AmoebaMinimizer &myMinimizer) {
    Vector startingPoint("TestPoint", numPC * 2);
    for (int i = 0; i < numPC * 2; ++i) {
        if (i < numPC)
            startingPoint[i] = PC[0][i];
        else
            startingPoint[i] = PC[1][i - numPC];
    }

    if(verbose) {
        std::cerr << "Start point:";
        for (int i = 0; i < numPC * 2; ++i) {
            std::cerr << startingPoint[i] << "\t";
        }
    }

    startingPoint.label = "startPoint";
    myMinimizer.func = &fn;
    myMinimizer.Reset(numPC * 2);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(epsilon);
    for (int i = 0; i < numPC; ++i) {
        PC[0][i] = myMinimizer.point[i];
    }
    for (int i = numPC; i < numPC * 2; ++i) {
        PC[1][i - numPC] = myMinimizer.point[i];
    }
    //fixAlpha usually converges well
    return true;
}

bool ContaminationEstimator::OptimizeHeterFixedPC(AmoebaMinimizer &myMinimizer) {
    return OptimizeHomo(myMinimizer);
}

bool ContaminationEstimator::OptimizeHomo(AmoebaMinimizer &myMinimizer) {
    Vector startingPoint("TestPoint", numPC + 1);
    for (int i = 0; i < numPC; ++i) {
        startingPoint[i] = PC[0][i];
    }
    startingPoint[numPC] = fullLLKFunc::Logit(alpha);
    if(verbose) {
        std::cerr << "Start point:";
        for (int i = 0; i < numPC; ++i) {
            std::cerr << startingPoint[i] << "\t";
        }
        std::cerr << "and alpha:\t" << alpha << std::endl;
    }
    startingPoint.label = "startPoint";
    myMinimizer.func = &fn;
    myMinimizer.Reset(numPC + 1);
    myMinimizer.point = startingPoint;
    double ret = myMinimizer.Minimize(epsilon);
    alpha = fullLLKFunc::invLogit(myMinimizer.point[numPC]);
    for (int i = 0; i < numPC; ++i) {
        PC[0][i] = myMinimizer.point[i];
    }
    if(ret == std::numeric_limits<double>::max()) return false;
    else return true;
}

bool ContaminationEstimator::OptimizeHomoFixedAlpha(AmoebaMinimizer &myMinimizer) {
    Vector startingPoint("TestPoint", numPC);
    for (int i = 0; i < numPC; ++i) {
        startingPoint[i] = PC[0][i];
    }
    if(verbose) {
        std::cerr << "Start point:";
        for (int i = 0; i < numPC; ++i) {
            std::cerr << startingPoint[i] << "\t";
        }
    }

    startingPoint.label = "startPoint";
    myMinimizer.func = &fn;
    myMinimizer.Reset(numPC);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(epsilon);
    for (int i = 0; i < numPC; ++i) {
        PC[0][i] = myMinimizer.point[i];
    }
    //fixAlpha usually converges well
    return true;
}

bool ContaminationEstimator::OptimizeHomoFixedPC(AmoebaMinimizer &myMinimizer) {
    Vector startingPoint("TestPoint", 1);
    startingPoint[0] = fullLLKFunc::Logit(alpha);

    if(verbose) {
        std::cerr << "Start point";
        std::cerr << "alpha:\t" << alpha << std::endl;
    }

    startingPoint.label = "startPoint";
    myMinimizer.func = &fn;
    myMinimizer.Reset(1);
    myMinimizer.point = startingPoint;
    double ret=myMinimizer.Minimize(epsilon);
    alpha = fullLLKFunc::invLogit(myMinimizer.point[0]);
    if(ret == std::numeric_limits<double>::max()) return false;//not converge
    else return true;
}

int ContaminationEstimator::ReadSVDMatrix(const std::string UDpath, const std::string Mean) {
    ReadMatrixUD(UDpath);
    ReadMean(Mean);
    return 0;
}

ContaminationEstimator::ContaminationEstimator(int nPC, const char *bedFile, int nThread, double ep)
        :
        numPC(nPC), PC(2, std::vector<PCtype>(nPC, 0.)), fn(nPC, this),numThread(nThread), epsilon(ep) {
    isAFknown = false;
    isPCFixed = false;
    isAlphaFixed = false;
    isHeter = true;
    ReadChooseBed(std::string(bedFile));
    alpha = 0.5;
    NumMarker = 0;

}

int ContaminationEstimator::ReadMatrixUD(const std::string &path) {
    std::ifstream fin(path);
    std::string line;
    std::vector<PCtype> tmpUD(numPC, 0);
    if (!fin.is_open()) {
        std::cerr << "Open file:" << path << "\t failed, exit!";
        exit(EXIT_FAILURE);
    }
    while (std::getline(fin, line)) {
        std::stringstream ss(line);
        //std::string chr;
        //int pos;
        //ss >> chr >> pos;
        for (int index = 0; index != numPC; ++index)
            ss >> tmpUD[index];
        UD.push_back(tmpUD);
        //initialize arrays
        NumMarker++;
        AFs.push_back(0.);
    }
    AF2s.assign(AFs.begin(), AFs.end());
    fin.close();
    return 0;
}

int ContaminationEstimator::ReadChooseBed(const std::string &path) {
    std::ifstream fin(path);
    std::string line, chr;
    int index(0), pos(0);
    char ref(0), alt(0);

    if (!fin.is_open()) {
        std::cerr << "Open file:" << path << "\t failed, exit!";
        exit(EXIT_FAILURE);
    }
    while (std::getline(fin, line)) {
        index++;
        std::stringstream ss(line);
        //std::string chr;
        //int pos;
        ss >> chr >> pos >> pos;
        ss >> ref >> alt;

        BedVec.push_back(region_t(chr, pos - 1, pos));
        PosVec.push_back(make_pair(chr, pos));
        ChooseBed[chr][pos] = std::make_pair(ref, alt);
    }
    fin.close();
    return 0;
}

int ContaminationEstimator::ReadMean(const std::string &path) {
    std::ifstream fin(path);
    std::string line;
    //int pos(0);
    double mu(0);
    std::string snpName, chr;
    if (!fin.is_open()) {
        std::cerr << "Open file:" << path << "\t failed, exit!";
        exit(EXIT_FAILURE);
    }
    while (std::getline(fin, line)) {
        std::stringstream ss(line);
        ss >> snpName;
        chr = snpName.substr(0, snpName.find(':', 0));
        //pos = atoi(snpName.substr(snpName.find(':', 0)+1,snpName.find('_',0)).c_str());
        ss >> mu;
        //std::cerr << chr << "\t" << pos << "\t"<<mu<<std::endl;
        //PosVec.push_back(make_pair(chr, pos));
        means.push_back(mu);
    }
    fin.close();
    return 0;
}

int ContaminationEstimator::ReadAF(const std::string &path) {
    std::ifstream fin(path);
    std::string line;
    uint32_t pos(0);
    double AF(0);
    std::string chr;
    char ref(0),alt(0);
//    int beg(0),end(0);
    if (!fin.is_open()) {
        std::cerr << "Open file:" << path << "\t failed, exit!";
        exit(EXIT_FAILURE);
    }
    while (std::getline(fin, line)) {
//        if(line[0]=='#'||line.find("INDEL")!=std::string::npos) continue;
        std::stringstream ss(line);
        ss >> chr;
        ss >> pos >>pos;
	ss >> ref >>alt;
//        ss>>snpName>>snpName>>snpName>>snpName>>snpName>>snpName;
//        beg=snpName.find("EUR_AF=");
//        beg+=7;
//        AF=atof(snpName.substr(beg,4).c_str());
        ss >> AF;
        knownAF[chr][pos] = AF;
    }
    return 0;
}

int ContaminationEstimator::ReadBam(const char *bamFile, const char *faiFile, const char *bedFile)
{
    viewer = SimplePileupViewer(&BedVec, bamFile, faiFile, bedFile, 1);
}