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
    double optimalPC1=0;
    double optimalPC2=0;
    double optimalPC3=0;
    double optimalPC4=0;
    double optimalAlpha =0;
    std::cout << "PCs in OptimizaLLK():" << std::endl;
    if(isPCFixed)
    {
        Vector startingPoint("TestPoint", 1);
        startingPoint[0] = alpha;
        startingPoint.label = "startPoint";
        myMinimizer.func = &fn;
        myMinimizer.Reset(1);
        myMinimizer.point = startingPoint;
        myMinimizer.Minimize(1e-6);
        optimalAlpha = fullLLKFunc::invLogit(myMinimizer.point[0]);
    }
    else if(isAlphaFixed)
    {
        Vector startingPoint("TestPoint", 2);
        startingPoint[0] = PC[0][0];
        startingPoint[1] = PC[0][1];
        startingPoint.label = "startPoint";
        myMinimizer.func = &fn;
        myMinimizer.Reset(2);
        myMinimizer.point = startingPoint;
        myMinimizer.Minimize(1e-6);
        optimalPC1 = myMinimizer.point[0];
        optimalPC2 = myMinimizer.point[1];
    }
    else {
        Vector startingPoint("TestPoint", 5);
        startingPoint[0] = PC[0][0];
        startingPoint[1] = PC[0][1];
        startingPoint[2] = PC[0][2];
        startingPoint[3] = PC[0][3];
        startingPoint[4] = alpha;
        startingPoint.label = "startPoint";
        myMinimizer.func = &fn;
        myMinimizer.Reset(5);
        myMinimizer.point = startingPoint;
        myMinimizer.Minimize(1e-6);
        optimalPC1 = myMinimizer.point[0];
        optimalPC2 = myMinimizer.point[1];
        optimalPC3 = myMinimizer.point[2];
        optimalPC4 = myMinimizer.point[3];
        optimalAlpha = fullLLKFunc::invLogit(myMinimizer.point[4]);
        std::cout << "PC3:" << optimalPC3 << "\tPC4:" << optimalPC4 << std::endl;
    }

    std::cout << "PC1:" << optimalPC1 << "\tPC2:" << optimalPC2 << std::endl;
    std::cout << "Alpha:" << (optimalAlpha<0.5?optimalAlpha:(1-optimalAlpha))<<std::endl;
    return 0;
}

int ContaminationEstimator::ReadSVDMatrix(const std::string UDpath, const std::string Mean, const std::string &Bed)
{
    ReadMatrixUD(UDpath);
    ReadMean(Mean);
    fn.initialize(this);
    return 0;
}

ContaminationEstimator::ContaminationEstimator(const char *bamFile, const char *faiFile, const char *bedFile) : PC(1, std::vector<PCtype>(4, 0)) {
    isAFknown = false;
    isPCFixed = false;
    isAlphaFixed = false;
    ReadChooseBed(std::string(bedFile));
    viewer = SimplePileupViewer(&BedVec,bamFile,faiFile,bedFile,1);
    alpha=0.01;
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
