/*
 * ContaminationEstimator.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: fanzhang
 */

#include "ContaminationEstimator.h"
#include <iostream>
#include <fstream>
#include <sstream>

#ifndef MPU_PATH
#define MPU_PATH "mpuTools"
#endif
ContaminationEstimator::ContaminationEstimator() {
	// TODO Auto-generated constructor stub

}

ContaminationEstimator::~ContaminationEstimator() {
	// TODO Auto-generated destructor stub
}


int ContaminationEstimator::OptimizeLLK()
{
    //std::cerr << "Now the label is:1" << std::endl;
    AmoebaMinimizer myMinimizer;
    //std::cerr << "Now the label is:2" << std::endl;
    Vector startingPoint("TestPoint",3);
    //std::cerr << "Now the label is:3" << std::endl;
    startingPoint[0] = PC[0][0];
    startingPoint[1] = PC[0][1];
    startingPoint[2] = alpha;
    //startingPoint.label = "startPoint";
    //std::cerr << "Now the label is:" << startingPoint.label << std::endl;
    myMinimizer.func = &fn;
    //std::cerr << "Before minimizing:1" << startingPoint.label << std::endl;
    myMinimizer.Reset(3);
    //std::cerr << "Before minimizing:2" << startingPoint.label << std::endl;
    myMinimizer.point = startingPoint;
    //std::cerr << "Before minimizing:3" << startingPoint.label << std::endl;
    myMinimizer.Minimize(1e-6);
    double optimalPC1 = myMinimizer.point[0];
    double optimalPC2 = myMinimizer.point[1];//fullLLKFunc::invLogit
    double optimalAlpha = fullLLKFunc::invLogit(myMinimizer.point[2]);
    std::cout << "PCs in OptimizaLLK():" << std::endl;
    std::cout << "PC1:" << optimalPC1 << "\tPC2:" << optimalPC2 << std::endl;
    std::cout << "Alpha:" << (optimalAlpha<0.5?optimalAlpha:(1-optimalAlpha))<<std::endl;

    return 0;
}
int ContaminationEstimator::RunFromVCF(const std::string VcfSiteAFFile,const std::string CurrentMPU, const std::string ReadGroup, const std::string Prefix)
{
	char cmdline[9*1024];
	char** params = new char*[9];
	sprintf(cmdline, MPU_PATH " verify --vcf %s --mpu %s --smID %s --out %s.ctm", VcfSiteAFFile.c_str(), CurrentMPU.c_str(), ReadGroup.c_str(), Prefix.c_str());
	std::stringstream ss(cmdline);
	std::string para;
	ss >> para;
	for (int i = 0; i != 9; ++i)
	{
		ss >> para;
		params[i] = new char[1024];
		strcpy(params[i], para.c_str());
	}
	//runVerify(9, params);
	return 0;
}

int ContaminationEstimator::ReadSVDMatrix(const std::string UDpath, const std::string Mean, const std::string &Bed)
{
    ReadMatrixUD(UDpath);
    ReadMean(Mean);
    fn.initialize(this);
    return 0;
}

ContaminationEstimator:: ContaminationEstimator(const char *bamFile, const char *faiFile, const char *bedFile, int nfiles):PC(1,std::vector<PCtype>(2,0)) {
    ReadChooseBed(std::string(bedFile));
    viewer = SimplePileupViewer(&BedVec,bamFile,faiFile,bedFile,1);
    alpha=0.01;
    NumMarker=0;
    NumIndividual=0;

}

int ContaminationEstimator::ReadMatrixUD(const std::string &path)
{
    std::ifstream fin(path);
    std::string line;
    uint32_t index(0);
    std::vector<PCtype> tmpUD(2, 0);
    if (!fin.is_open()) {  std::cerr<<"Open file:"<<path<<"\t failed, exit!";exit(EXIT_FAILURE);  }
    while (std::getline(fin, line))
    {
        std::stringstream ss(line);
        //std::string chr;
        //int pos;
        //ss >> chr >> pos;
        ss >> tmpUD[0] >> tmpUD[1];
        UD.push_back(tmpUD);
        //initialize arrays
        NumMarker++;
        AFs.push_back(0.);
    }
    fin.close();
    return UD.size();
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
        ChooseBed[chr][pos] = std::make_pair(ref,alt);
    }
    fin.close();
    return UD.size();
}

int ContaminationEstimator::ReadMean(const std::string &path)
{
    std::ifstream fin(path);
    std::string line;
    int index(0),pos(0);
    double mu(0);
    std::string snpName,chr;
    if (!fin.is_open()) {  std::cerr<<"Open file:"<<path<<"\t failed, exit!";exit(EXIT_FAILURE);  }
    while (std::getline(fin, line))
    {
        std::stringstream ss(line);
        ss >> snpName;
        chr=snpName.substr(0,snpName.find(':',0));
        pos = atoi(snpName.substr(snpName.find(':', 0)+1,snpName.find('_',0)).c_str());
        ss >> mu;
        //std::cerr << chr << "\t" << pos << "\t"<<mu<<std::endl;
        PosVec.push_back(make_pair(chr, pos));
        means.push_back(mu);
        //means[index]=mu;
        //index++;
    }
    fin.close();
    return means.size();
}