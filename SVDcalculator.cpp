#include <unordered_map>
#include <fstream>
#include <Error.h>
#include <unordered_set>
#include "SVDcalculator.h"
#include "libVcf/libVcfVcfFile.h"
#include "Eigen/Dense"
using namespace libVcf;
using namespace Eigen;

SVDcalculator::SVDcalculator()
{
    numIndividual = 0;
    numMarker = 0;
}

SVDcalculator::~SVDcalculator() {}

int SVDcalculator::ReadVcf(const std::string &VcfPath,
                           std::vector<std::vector<char> >& genotype,
                           int & nSamples, int& nMarkers) {
    try {
        int maxPhred=255;
        VcfFile *pVcf = new VcfFile;
        pVcf->bSiteOnly = false;
        pVcf->bParseGenotypes = false;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = true;
        pVcf->openForRead(VcfPath.c_str());

        // check the sanity of data
        if (pVcf->getSampleCount() == 0) {
            throw VcfFileException("No individual genotype information exist in the input VCF file %s",
                                   VcfPath.c_str());
        }

        std::unordered_set<std::string> acceptChr={"1","2","3","4","5","6","7","8","9","10",
                                              "11","12","13","14","15","16","17","18","19","20",
                                              "21","22",
                                              "chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chr8", "chr9","chr10",
                                              "chr11","chr12","chr13","chr14","chr15","chr16","chr17", "chr18","chr19",
                                              "chr20","chr21","chr22"};

        nSamples = pVcf->getSampleCount();
        nMarkers = 0;
        char refAllele;
        char altAllele;
        VcfMarker *pMarker = new VcfMarker;
        String markerName;
        String prevMarkerName;
        while (pVcf->iterateMarker()) {//for each marker

            pMarker = pVcf->getLastMarker();
            markerName.printf("%s:%d", pMarker->sChrom.c_str(), pMarker->nPos);
            if(prevMarkerName==markerName)
            {
                fprintf(stderr,"Duplicated Marker at %s\n",markerName.c_str());
                exit(EXIT_FAILURE);
            }
            if(acceptChr.find(std::string(pMarker->sChrom.c_str())) == acceptChr.end())
            {
                fprintf(stderr,"skip non-autosome at %s\n",markerName.c_str());
                continue;
            }
//            if(pMarker->sChrom.Compare("X")==0 or pMarker->sChrom.Compare("chrX")==0 or
//              pMarker->sChrom.Compare("Y")==0 or pMarker->sChrom.Compare("chrY")==0)
//            {
//                fprintf(stderr,"skip non-autosome at %s\n",markerName.c_str());
//                continue;
//            }
            if(pMarker->sRef.Length()>1 or pMarker->asAlts[0].Length()>1 )
            {
                fprintf(stderr,"skip indel at %s\n",markerName.c_str());
                continue;
            }
            refAllele=pMarker->sRef[0];
            altAllele=pMarker->asAlts[0][0];

            chooseBed[pMarker->sChrom.c_str()][pMarker->nPos]=std::make_pair(refAllele,altAllele);
            BedVec.push_back(region_t(pMarker->sChrom.c_str(),pMarker->nPos-1,pMarker->nPos));
            int PLidx = pMarker->asFormatKeys.Find("PL");
            int PLGLGTflag = 0;//0 for PL, 1 for GL, 2 for GT
            if (PLidx < 0) {
                PLidx = pMarker->asFormatKeys.Find("GL");
                if (PLidx >= 0) PLGLGTflag = 1;//found GL
                else {
                    PLidx = pMarker->asFormatKeys.Find("GT");
                    if(PLidx >= 0) PLGLGTflag = 2;//found GT
                    else throw VcfFileException("Cannot recognize GT, GL or PL key in FORMAT field");
                }
            }
            //printf("reading vcf 1\n\n");
            int formatLength = pMarker->asFormatKeys.Length();
            int idx11 = 0, idx12 = 1, idx22 = 2;

            StringArray phred;

            long phred11;
            long phred12;
            long phred22;
            std::vector<char> perMarkerGeno(nSamples,-1);
            for (int i = 0; i < nSamples; i++)//for each individual
            {
                    if(nMarkers==0) Samples.push_back(pVcf->getSampleID(i).c_str());
                    if(PLGLGTflag==0)//found PL
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i * formatLength], ",");
                        phred11=phred[idx11].AsInteger();
                        phred12=phred[idx12].AsInteger();
                        phred22=phred[idx22].AsInteger();
                    }
                    else if(PLGLGTflag == 1)//found GL
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i * formatLength], ",");
                        phred11=static_cast<int>(-10. * phred[idx11].AsDouble());
                        phred12=static_cast<int>(-10. * phred[idx12].AsDouble());
                        phred22=static_cast<int>(-10. * phred[idx22].AsDouble());
                    }
                    else//found GT
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i * formatLength], "|/");
                        long geno=phred[0].AsInteger()+phred[1].AsInteger();
                        if(geno==0)
                        {
                            phred11=0;
                            phred12=30;
                            phred22=50;
                        }
                        else if(geno==1)
                        {
                            phred11=50;
                            phred12=0;
                            phred22=50;
                        }
                        else
                        {
                            phred11=50;
                            phred12=30;
                            phred22=0;
                        }
                    }

                    if ((phred11 < 0) || (phred12 < 0) || (phred22 < 0)) {
                        error("Negative PL or Positive GL observed");
                    }

                    if (phred11 > maxPhred) phred11 = maxPhred;
                    if (phred12 > maxPhred) phred12 = maxPhred;
                    if (phred22 > maxPhred) phred22 = maxPhred;
//
//                    printf("phred scores are %f, %f, %f;\tphred11/12/22 %d, %d, %d\n", phred[idx11].AsDouble(), phred[idx12].AsDouble(), phred[idx22].AsDouble(),phred11,phred12,phred22);
                    int minGeno = -1;
                    int minPhred = maxPhred;
                    if(phred11 < minPhred)
                    {
                        minPhred = phred11;
                        minGeno = 0;
                    }
                    if(phred12 < minPhred)
                    {
                        minPhred = phred12;
                        minGeno = 1;
                    }
                    if(phred22 < minPhred)
                    {
                        minPhred = phred22;
                        minGeno = 2;
                    }
                    perMarkerGeno[i] = minGeno;

//                    fprintf(stderr,"marker:%d\t%d\t%d\t%d\n",markerindex,engine.genotypes[personIndices[i]][genoindex],engine.genotypes[personIndices[i]][genoindex + 1],engine.genotypes[personIndices[i]][genoindex + 2] );
            }
            genotype.push_back(perMarkerGeno);
            nMarkers++;
            prevMarkerName = markerName;
        }

        delete pVcf;
        //delete pMarker;
    }
    catch (VcfFileException& e) {
        error(e.what());
    }
    return 0;
}

void SVDcalculator::ProcessRefVCF(const std::string &VcfPath)
{

    std::vector<std::vector<char> > genotype;//markers X samples

    ReadVcf(VcfPath, genotype, numIndividual, numMarker);
    MatrixXf genoMatrix(numMarker,numIndividual);
//    std::cerr<<numMarker<<"\t"<<genotype.size()<<"\t"<<numIndividual<<"\t"<<genotype[0].size()<<std::endl;
    notice("Number of Markers:%d\n",numMarker);
    notice("Number of Individuals:%d\n",numIndividual);

    for (int i = 0; i <genotype.size() ; ++i) {//per marker
        for (int j = 0; j <genotype[i].size() ; ++j) {//per sample
            genoMatrix(i,j)=genotype[i][j];
        }
    }
    genotype.clear();
    VectorXf mu = genoMatrix.rowwise().mean();
    for (int rowIdx = 0; rowIdx < numMarker; ++rowIdx) {
        for (int colIdx = 0; colIdx < numIndividual; ++colIdx) {
            genoMatrix(rowIdx,colIdx)-=mu(rowIdx);
        }
        Mu.push_back(mu(rowIdx));
    }
    JacobiSVD<MatrixXf> svd(genoMatrix, ComputeThinU | ComputeThinV);
    auto matrixD = svd.singularValues().asDiagonal();
    MatrixXf matrixUD = svd.matrixU() * matrixD;//marker X PC
    UD.resize(matrixUD.rows(),std::vector<double>(matrixUD.cols(),0.f));
    for (int rowIdx = 0; rowIdx <matrixUD.rows(); ++rowIdx) {
        for (int colIdx = 0; colIdx <matrixUD.cols(); ++colIdx) {
            UD[rowIdx][colIdx] = matrixUD(rowIdx,colIdx);
        }
    }
    MatrixXf matrixV = svd.matrixV();
    PC.resize(numIndividual,std::vector<double>(matrixUD.cols(),0.f));
    for (int sampleIdx = 0; sampleIdx < numIndividual ; ++sampleIdx) {
        for (int pcIdx = 0; pcIdx <matrixUD.cols() ; ++pcIdx) {
            PC[sampleIdx][pcIdx] = matrixV(sampleIdx,pcIdx);
        }
    }

    WriteSVD(VcfPath);
}

vector<vector<double>> SVDcalculator::GetUDMatrix() {
    return UD;
}

vector<vector<double>> SVDcalculator::GetPCMatrix() {
    return PC;
}

std::vector<PCtype> SVDcalculator::GetMuArray() {
    return Mu;
}

BED SVDcalculator::GetchooseBed() {
    return chooseBed;
}

vector<region_t> SVDcalculator::GetBedVec() {
    return BedVec;
}

void SVDcalculator::WriteSVD(const std::string &Prefix) {
    std::ofstream fMu(Prefix+".mu");
    std::ofstream fUD(Prefix+".UD");
    std::ofstream fPC(Prefix+".PC");
    std::ofstream fBed(Prefix+".bed");
    std::string chr;
    int beg(0),end(0);
    for (int i = 0; i < numMarker; ++i) {
        chr=BedVec[i].chr;
        beg=BedVec[i].beg;
        end=BedVec[i].end;
        fMu<<chr+":"+std::to_string(end)<<"\t"<<Mu[i]<<std::endl;
        fBed<<chr<<"\t"<<beg<<"\t"<<end<<"\t"<<chooseBed[chr][end].first<<"\t"<<chooseBed[chr][end].second<<std::endl;
        for (int j = 0; j < UD[i].size() ; ++j) {
            fUD<<UD[i][j]<<"\t";
        }
        fUD<<std::endl;
    }
    for (int k = 0; k <numIndividual; ++k) {
        fPC<<Samples[k]<<"\t";
        for (int i = 0; i < PC[k].size(); ++i) {
            fPC<<PC[k][i]<<"\t";
        }
        fPC<<std::endl;
    }
    fBed.close();
    fMu.close();
    fUD.close();
    fPC.close();
}
