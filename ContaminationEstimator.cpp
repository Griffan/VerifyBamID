/*
 * ContaminationEstimator.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: fanzhang
 */

#include "ContaminationEstimator.h"
#include <fstream>
#include <sstream>

#ifndef MPU_PATH
#define MPU_PATH "mpuTools"
#endif
ContaminationEstimator::ContaminationEstimator() {

}

ContaminationEstimator::~ContaminationEstimator() {
}


int ContaminationEstimator::OptimizeLLK()
{
    AmoebaMinimizer myMinimizer;

    if(!isHeter)
    {
        if(isPCFixed) {
            std::cout << "Estimation from OptimizeHomFixedPC:"<<std::endl;
            OptimizeHomFixedPC(myMinimizer);
        }
        else if(isAlphaFixed) {
            std::cout << "Estimation from OptimizeHomFixedAlpha:"<<std::endl;
            OptimizeHomFixedAlpha(myMinimizer);
        }
        else
        {
            for (int k = 0; k <PC[0].size(); ++k) {
                PC[0][k]=static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            }
            alpha = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);

            std::cout << "Estimation from OptimizeHom:"<<std::endl;
            OptimizeHom(myMinimizer);
        }
        std::cout << "PC1:" << PC[0][0] << "\tPC2:" << PC[0][1] << std::endl;
    }
    else//contamination source from different population
    {
        if(isPCFixed) {
            std::cout << "Estimation from OptimizeHeterFixedPC:"<<std::endl;
            OptimizeHeterFixedPC(myMinimizer);
        }
        else if(isAlphaFixed) {
            std::cout << "Estimation from OptimizeHeterFixedAlpha:"<<std::endl;
            OptimizeHeterFixedAlpha(myMinimizer);
        }
        else
        {
            for (int k = 0; k <PC[0].size(); ++k) {
                PC[0][k]=static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            }
            for (int k = 0; k <PC[1].size(); ++k) {
                PC[1][k]=static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            }
            alpha = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            std::cout << "Estimation from OptimizeHeter:"<<std::endl;
            isHeter=false;
            OptimizeHom(myMinimizer);
            std::cerr << "testAlpha:" << (alpha<0.5?alpha:(1-alpha))<<std::endl;
            isHeter=true;
            fn.initialize();
            OptimizeHeter(myMinimizer);
        }
        std::cout << "PC1:" << PC[0][0] << "\tPC2:" << PC[0][1] << std::endl;
        std::cout << "PC3:" << PC[1][0] << "\tPC4:" << PC[1][1] << std::endl;
    }

    std::cout << "Alpha:" << (alpha<0.5?alpha:(1-alpha))<<std::endl;
    return 0;
}

void ContaminationEstimator::OptimizeHeter(AmoebaMinimizer &myMinimizer) {
    Vector startingPoint("TestPoint", numPC * 2 + 1);
    for (int i = 0; i < numPC * 2; ++i) {
                if(i < numPC)
                    startingPoint[i] = PC[0][i];
                else
                    startingPoint[i] = PC[1][i-numPC];
            }
    startingPoint[numPC * 2] = alpha;
    startingPoint.label = "startPoint";
    myMinimizer.func = &fn;
    myMinimizer.Reset(numPC * 2 + 1);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(1e-6);
    alpha = fullLLKFunc::invLogit(myMinimizer.point[numPC * 2]);
    for (int i = 0; i < numPC; ++i) {
        PC[0][i] = myMinimizer.point[i];
    }
    for (int i = numPC; i < numPC * 2; ++i) {
        PC[1][i-numPC] = myMinimizer.point[i];
    }
}

void ContaminationEstimator::OptimizeHeterFixedAlpha(AmoebaMinimizer &myMinimizer) {
    Vector startingPoint("TestPoint", numPC * 2);
    for (int i = 0; i < numPC * 2; ++i) {
                if(i < numPC)
                    startingPoint[i] = PC[0][i];
                else
                    startingPoint[i] = PC[1][i-numPC];
            }
    startingPoint.label = "startPoint";
    myMinimizer.func = &fn;
    myMinimizer.Reset(numPC * 2);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(1e-6);
    for (int i = 0; i < numPC; ++i) {
        PC[0][i] = myMinimizer.point[i];
    }
    for (int i = numPC; i < numPC * 2; ++i) {
        PC[1][i-numPC] = myMinimizer.point[i];
    }
}

void ContaminationEstimator::OptimizeHeterFixedPC(AmoebaMinimizer &myMinimizer) {
    OptimizeHomFixedPC(myMinimizer);
}

void ContaminationEstimator::OptimizeHom(AmoebaMinimizer &myMinimizer) {
    Vector startingPoint("TestPoint", numPC + 1);
    for (int i = 0; i < numPC; ++i) {
                startingPoint[i] = PC[0][i];
            }
    startingPoint[numPC] = alpha;
    startingPoint.label = "startPoint";
    myMinimizer.func = &fn;
    myMinimizer.Reset(numPC + 1);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(1e-6);
    alpha = fullLLKFunc::invLogit(myMinimizer.point[numPC]);
    for (int i = 0; i < numPC; ++i) {
        PC[0][i] = myMinimizer.point[i];
    }
}

void ContaminationEstimator::OptimizeHomFixedAlpha(AmoebaMinimizer &myMinimizer) {
    Vector startingPoint("TestPoint", numPC);
    for (int i = 0; i < numPC; ++i) {
                startingPoint[i] = PC[0][i];
            }
    startingPoint.label = "startPoint";
    myMinimizer.func = &fn;
    myMinimizer.Reset(2);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(1e-6);
    for (int i = 0; i < numPC; ++i) {
        PC[0][i] = myMinimizer.point[i];
    }
}

void ContaminationEstimator::OptimizeHomFixedPC(AmoebaMinimizer &myMinimizer) {
    Vector startingPoint("TestPoint", 1);
    startingPoint[0] = alpha;
    startingPoint.label = "startPoint";
    myMinimizer.func = &fn;
    myMinimizer.Reset(1);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(1e-6);
    alpha = fullLLKFunc::invLogit(myMinimizer.point[0]);
}

int ContaminationEstimator::ReadSVDMatrix(const std::string UDpath, const std::string Mean, const std::string &Bed)
{
    ReadMatrixUD(UDpath);
    ReadMean(Mean);
    fn.initialize();
    return 0;
}

ContaminationEstimator::ContaminationEstimator(const char *bamFile, const char *faiFile, const char *bedFile) :
        numPC(2), PC(2, std::vector<PCtype>(numPC, 0)),fn(numPC,this) {
    isAFknown = false;
    isPCFixed = false;
    isAlphaFixed = false;
    isHeter = true;

    ReadChooseBed(std::string(bedFile));
    viewer = SimplePileupViewer(&BedVec,bamFile,faiFile,bedFile,1);
    alpha=0.99;
    NumMarker=0;

}

int ContaminationEstimator::ReadMatrixUD(const std::string &path)
{
    std::ifstream fin(path);
    std::string line;
    std::vector<PCtype> tmpUD(4, 0);
    if (!fin.is_open()) {  std::cerr<<"Open file:"<<path<<"\t failed, exit!";exit(EXIT_FAILURE);  }
    while (std::getline(fin, line))
    {
        std::stringstream ss(line);
        //std::string chr;
        //int pos;
        //ss >> chr >> pos;
        ss >> tmpUD[0] >> tmpUD[1]>>tmpUD[2]>>tmpUD[3];
        UD.push_back(tmpUD);
        //initialize arrays
        NumMarker++;
        AFs.push_back(0.);
    }
    AF2s.assign(AFs.begin(),AFs.end());
    fin.close();
    return 0;
}

int ContaminationEstimator::ReadChooseBed(const std::string &path)
{
    std::ifstream fin(path);
    std::string line,chr;
    uint32_t index(0),pos(0);
    char ref(0), alt(0);

    if (!fin.is_open()) {  std::cerr<<"Open file:"<<path<<"\t failed, exit!";exit(EXIT_FAILURE);  }
    while (std::getline(fin, line))
    {
        index++;
        std::stringstream ss(line);
        //std::string chr;
        //int pos;
        ss >> chr >> pos>>pos;
        ss >> ref >> alt;

        BedVec.push_back(region_t(chr,pos-1,pos));
        PosVec.push_back(make_pair(chr, pos));
	ChooseBed[chr][pos] = std::make_pair(ref,alt);
    }
    fin.close();
    return 0;
}

int ContaminationEstimator::ReadMean(const std::string &path)
{
    std::ifstream fin(path);
    std::string line;
    //int pos(0);
    double mu(0);
    std::string snpName,chr;
    if (!fin.is_open()) {  std::cerr<<"Open file:"<<path<<"\t failed, exit!";exit(EXIT_FAILURE);  }
    while (std::getline(fin, line))
    {
        std::stringstream ss(line);
        ss >> snpName;
        chr=snpName.substr(0,snpName.find(':',0));
        //pos = atoi(snpName.substr(snpName.find(':', 0)+1,snpName.find('_',0)).c_str());
        ss >> mu;
        //std::cerr << chr << "\t" << pos << "\t"<<mu<<std::endl;
        //PosVec.push_back(make_pair(chr, pos));
        means.push_back(mu);
    }
    fin.close();
    return 0;
}

int ContaminationEstimator::ReadAF(const std::string & path)
{
    std::ifstream fin(path);
    std::string line;
    uint32_t pos(0);
    double AF(0);
    std::string chr;
//    int beg(0),end(0);
    if (!fin.is_open()) {  std::cerr<<"Open file:"<<path<<"\t failed, exit!";exit(EXIT_FAILURE);  }
    while (std::getline(fin, line))
    {
//        if(line[0]=='#'||line.find("INDEL")!=std::string::npos) continue;
        std::stringstream ss(line);
        ss>>chr;
        ss>>pos;
//        ss>>snpName>>snpName>>snpName>>snpName>>snpName>>snpName;
//        beg=snpName.find("EUR_AF=");
//        beg+=7;
//        AF=atof(snpName.substr(beg,4).c_str());
        ss>>AF;
        knownAF[chr][pos]=AF;
    }
    return 0;
}
