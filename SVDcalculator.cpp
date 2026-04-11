#include "SVDcalculator.h"
#include "Eigen/Dense"
#include "libVcf/libVcfFile.h"
#include <Error.h>
#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
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
                           int & nSamples, int& nMarkers,
                           const std::unordered_set<std::string>& includeChr) {
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

        bool filterByChrom = !includeChr.empty();
        if (filterByChrom) {
            notice("Filtering to %d chromosome(s) specified by --IncludeChr", (int)includeChr.size());
        }

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
                error("Duplicated Marker: %s",markerName.c_str());
            }
            if(pMarker->asFilters.Length()>1 || pMarker->asFilters[0]!="PASS")
            {
              warning("Skip filtered (%s) marker: %s",pMarker->asFilters[0].c_str(), markerName.c_str());
              continue;
            }
            if(pMarker->asAlts.Length()>1)
            {
              warning("Skip non-Biallelic marker: %s",markerName.c_str());
              continue;
            }
            if(pMarker->sRef.Length()>1 or pMarker->asAlts[0].Length()>1 )
            {
                warning("Skip non-SNP marker: %s",markerName.c_str());
                continue;
            }
            if(filterByChrom && includeChr.find(std::string(pMarker->sChrom.c_str())) == includeChr.end())
            {
              continue;
            }
            refAllele=pMarker->sRef[0];
            altAllele=pMarker->asAlts[0][0];

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
            int formatLength = pMarker->asFormatKeys.Length();
            int idx11 = 0, idx12 = 1, idx22 = 2;

            StringArray phred;

            long phred11 = 0;
            long phred12 = 0;
            long phred22 = 0;
            std::vector<char> perMarkerGeno(nSamples,-1);
            int nMissingGenoSamples = 0;
            for (int i = 0; i < nSamples; i++)//for each individual
            {
                    if(nMarkers==0) Samples.push_back(pVcf->getSampleID(i).c_str());
                    if(PLGLGTflag==0)//found PL
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i * formatLength], ",");
                        if (phred.Length() == 3 && phred[0] != ".") {
                          phred11 = phred[idx11].AsInteger();
                          phred12 = phred[idx12].AsInteger();
                          phred22 = phred[idx22].AsInteger();
                        }
                        else nMissingGenoSamples++;
                    }
                    else if(PLGLGTflag == 1)//found GL
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i * formatLength], ",");
                        if (phred.Length() == 3 && phred[0] != ".") {
                          phred11 =
                              static_cast<int>(-10. * phred[idx11].AsDouble());
                          phred12 =
                              static_cast<int>(-10. * phred[idx12].AsDouble());
                          phred22 =
                              static_cast<int>(-10. * phred[idx22].AsDouble());
                        }
                        else nMissingGenoSamples++;
                    }
                    else//found GT
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i * formatLength], "|/");
                        if (phred.Length() == 3 && phred[0] != ".") {
                          long geno =
                              phred[0].AsInteger() + phred[1].AsInteger();
                          if (geno == 0) {
                            phred11 = 0;
                            phred12 = 30;
                            phred22 = 50;
                          } else if (geno == 1) {
                            phred11 = 50;
                            phred12 = 0;
                            phred22 = 50;
                          } else {
                            phred11 = 50;
                            phred12 = 30;
                            phred22 = 0;
                          }
                        }
                        else nMissingGenoSamples++;
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
                    long minPhred = maxPhred;
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
            }

            float genoMissingRate = static_cast<float>(nMissingGenoSamples) / nSamples;
            if(genoMissingRate > 0.2f)
            {
              warning("Skip marker (%s) with high missing rate (%f > 0.2) in genotype fields.",markerName.c_str(), genoMissingRate);
              continue;
            }
            genotype.push_back(perMarkerGeno);
            chooseBed[pMarker->sChrom.c_str()][pMarker->nPos]=std::make_pair(refAllele,altAllele);
            BedVec.push_back(region_t(pMarker->sChrom.c_str(),pMarker->nPos-1,pMarker->nPos));

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

void SVDcalculator::ProcessRefVCF(const std::string &VcfPath,
                                  const std::unordered_set<std::string>& includeChr, 
                                  bool skipMinSampleCountCheck,
                                  int numSVDPCs)
{

    std::vector<std::vector<char> > genotype;//markers X samples

    ReadVcf(VcfPath, genotype, numIndividual, numMarker, includeChr);
    MatrixXf genoMatrix(numMarker,numIndividual);
    notice("Number of Markers:%d",numMarker);
    notice("Number of Individuals:%d",numIndividual);

    if(numMarker < 5000)
    {
      error("Insufficient number of markers (need >= 5000, have %d)\n", numMarker);
    }

    if(numIndividual < 1000)
    {
      if (skipMinSampleCountCheck) {
        warning("Only %d individuals in reference panel (recommended minimum is 1000). "
                "Proceeding because --SkipMinSampleCountCheck is set. Contamination "
                "estimates may be unreliable if the panel does not adequately capture "
                "population structure.", numIndividual);
      } else {
        error("Insufficient number of individuals (need >= 1000, have %d). If your "
              "reference panel adequately captures population structure with fewer "
              "samples, rerun with --SkipMinSampleCountCheck.\n", numIndividual);
      }
    }

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

    // Log proportion of variance explained by each principal component.
    // Variance explained by PC_i = sigma_i^2 / sum(sigma_j^2) where sigma
    // are the singular values. This helps users choose an appropriate --NumPC.
    auto singularValues = svd.singularValues();
    double totalVariance = 0.0;
    for (int i = 0; i < singularValues.size(); ++i) {
        totalVariance += (double)singularValues[i] * (double)singularValues[i];
    }
    notice("SVD completed. Total singular values: %d", (int)singularValues.size());
    double cumulative = 0.0;
    int numToLog = std::min((int)singularValues.size(), 20);
    if (totalVariance <= 0.0) {
        warning("Total variance is zero after SVD; skipping variance-explained logging for principal components.");
    } else {
        for (int i = 0; i < numToLog; ++i) {
            double sv = singularValues[i];
            double varExplained = (sv * sv) / totalVariance;
            cumulative += varExplained;
            notice("  PC%d: singular_value=%.4f  variance_explained=%.4f (%.2f%%)  cumulative=%.4f (%.2f%%)",
                   i + 1, sv, varExplained, varExplained * 100.0, cumulative, cumulative * 100.0);
        }
        if ((int)singularValues.size() > numToLog) {
            notice("  ... (%d more components not shown)", (int)singularValues.size() - numToLog);
        }
    }

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

    WriteSVD(VcfPath, numSVDPCs);
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

void SVDcalculator::WriteSVD(const std::string &Prefix, int numSVDPCs) {
    // Determine how many PCs to write
    int numAvailable = (numMarker > 0 && !UD.empty()) ? (int)UD[0].size() : 0;
    int numToWrite;
    if (numSVDPCs <= 0) {
        numToWrite = numAvailable;
    } else {
        numToWrite = std::min(numSVDPCs, numAvailable);
    }
    notice("Writing SVD output files with %d PCs (of %d available) to prefix: %s",
           numToWrite, numAvailable, Prefix.c_str());

    std::ofstream fMu(Prefix+".mu");
    std::ofstream fUD(Prefix+".UD");
    std::ofstream fPC(Prefix+".V");
    std::ofstream fBed(Prefix+".bed");
    std::string chr;
    int beg(0),end(0);
    for (int i = 0; i < numMarker; ++i) {
        chr=BedVec[i].chr;
        beg=BedVec[i].beg;
        end=BedVec[i].end;
        fMu<<chr+":"+std::to_string(end)<<"\t"<<Mu[i]<<std::endl;
        fBed<<chr<<"\t"<<beg<<"\t"<<end<<"\t"<<chooseBed[chr][end].first<<"\t"<<chooseBed[chr][end].second<<std::endl;
        for (int j = 0; j < numToWrite; ++j) {
            fUD<<UD[i][j]<<"\t";
        }
        fUD<<std::endl;
    }
    for (int k = 0; k <numIndividual; ++k) {
        fPC<<Samples[k]<<"\t";
        for (int i = 0; i < numToWrite; ++i) {
            fPC<<PC[k][i]<<"\t";
        }
        fPC<<std::endl;
    }
    fBed.close();
    fMu.close();
    fUD.close();
    fPC.close();
}
