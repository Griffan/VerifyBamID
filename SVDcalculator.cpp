//
// Created by Fan Zhang on 11/2/16.
//

#include <unordered_map>
#include "SVDcalculator.h"
#include "libVcf/libVcfVcfFile.h"
using namespace libVcf;
SVDcalculator::SVDcalculator()
{
    numIndividual = 0;
    numMarker = 0;
}

SVDcalculator::~SVDcalculator() {}

int SVDcalculator::ReadVcf(const std::string &VcfPath) {
    printf("starting LoadGenotypeFromUnphasedVCF\n\n");
    try {
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

        int nSamples = pVcf->getSampleCount();
        std::unordered_map<int, int> personIndices;


        int markerindex = 0;
        VcfMarker *pMarker = new VcfMarker;
        String markerName;
        while (pVcf->iterateMarker()) {//for each marker

            pMarker = pVcf->getLastMarker();
            markerName.printf("%s:%d", pMarker->sChrom.c_str(), pMarker->nPos);

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
            int genoindex = markerindex * 3;

            long phred11;
            long phred12;
            long phred22;
            for (int i = 0; i < nSamples; i++)//for each individual
            {
                if (personIndices.find(i) != personIndices.end()) {

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

                    engine.genotypes[personIndices[i]][genoindex] = phred11;
                    engine.genotypes[personIndices[i]][genoindex + 1] = phred12;
                    engine.genotypes[personIndices[i]][genoindex + 2] = phred22;
//                    fprintf(stderr,"marker:%d\t%d\t%d\t%d\n",markerindex,engine.genotypes[personIndices[i]][genoindex],engine.genotypes[personIndices[i]][genoindex + 1],engine.genotypes[personIndices[i]][genoindex + 2] );
                }
            }
        }

        delete pVcf;
        //delete pMarker;
    }
    catch (VcfFileException e) {
        error(e.what());
    }
    return 0;
}
