#include "ContaminationEstimator.h"
#include <chrono>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

// Simple RAII timer that logs elapsed wall-clock time for a named phase.
namespace {
struct PhaseTimer {
    std::string name;
    std::chrono::steady_clock::time_point start;

    explicit PhaseTimer(const std::string &phaseName)
        : name(phaseName), start(std::chrono::steady_clock::now()) {
        notice("  Starting phase: %s", name.c_str());
    }

    ~PhaseTimer() {
        auto end = std::chrono::steady_clock::now();
        double secs = std::chrono::duration<double>(end - start).count();
        notice("  Finished phase: %s  [%.3f seconds]", name.c_str(), secs);
    }
};
} // anonymous namespace

ContaminationEstimator::ContaminationEstimator() {

}

ContaminationEstimator::~ContaminationEstimator() {
}

ContaminationEstimator::ContaminationEstimator(int nPC, const char *bedFile, int nThread, double ep)
        :
        numPC(nPC), PC(2, std::vector<PCtype>(nPC, 0.)),muv(numPC,0), sdv(numPC,0), fn(nPC, this), numThread(nThread), epsilon(ep) {
    isAFknown = false;
    isPCFixed = false;
    isAlphaFixed = false;
    isHeter = true;
    isPileupInput = false;
    isSanityCheckDisabled = false;
    ReadChooseBed(std::string(bedFile));
    alpha = 0.5;
    NumMarker = 0;

}

// Resolve per-marker data that is constant across all optimization iterations.
//
// During optimization, ComputeMixLLKs is called tens of thousands of times.
// Each call iterates over all markers and previously performed nested hash map
// lookups (posIndex[chr][pos], ChooseBed[chr][pos], knownAF[chr][pos]) on every
// iteration.  Since none of these values change between calls, we resolve them
// once here into a flat vector indexed by marker ordinal.
//
// For markers not present in the viewer (e.g. the BAM had no reads at that
// position), baseInfoIndex is set to -1 so ComputeMixLLKs can skip them with
// a single integer comparison instead of two hash lookups.
//
// Must be called after data loading (ReadBam/ReadPileup) and before the
// optimization loop in OptimizeLLK.
void ContaminationEstimator::BuildResolvedMarkers() {
    resolvedMarkers.resize(NumMarker);
    for (size_t i = 0; i < NumMarker; ++i) {
        const std::string& chr = PosVec[i].first;
        int pos = PosVec[i].second;
        ResolvedMarker& rm = resolvedMarkers[i];

        auto chrIt = viewer.posIndex.find(chr);
        if (chrIt == viewer.posIndex.end()) { rm.baseInfoIndex = -1; continue; }
        auto posIt = chrIt->second.find(pos);
        if (posIt == chrIt->second.end()) { rm.baseInfoIndex = -1; continue; }

        rm.baseInfoIndex = posIt->second;
        rm.altBase = ChooseBed[chr][pos].second;
        rm.knownAFValue = 0.0;
        if (isAFknown) {
            rm.knownAFValue = knownAF[chr][pos];
        }
    }
}

int ContaminationEstimator::OptimizeLLK(const std::string &OutputPrefix) {
    AmoebaMinimizer myMinimizer;

    BuildResolvedMarkers();

    {
        PhaseTimer t("Initialize likelihood");
        fn.Initialize();
    }

    if (!isHeter) {
        if (isPCFixed) {
            std::cout << "Estimation from OptimizeHomoFixedPC:" << std::endl;
            PhaseTimer t("OptimizeHomoFixedPC");
            OptimizeHomoFixedPC(myMinimizer);
        } else if (isAlphaFixed) {
            PhaseTimer t("OptimizeHomoFixedAlpha");
            OptimizeHomoFixedAlpha(myMinimizer);
        } else {
            std::cout << "Estimation from OptimizeHomo:" << std::endl;
            PhaseTimer t("OptimizeHomo");
            OptimizeHomo(myMinimizer);
        }
    } else//contamination source from different population
    {
        if (isPCFixed) {
            std::cout << "Estimation from OptimizeHeterFixedPC:" << std::endl;
            PhaseTimer t("OptimizeHeterFixedPC");
            OptimizeHeterFixedPC(myMinimizer);
        } else if (isAlphaFixed) {
            std::cout << "Estimation from OptimizeHeterFixedAlpha:" << std::endl;
            {
                PhaseTimer t("OptimizeHomoFixedAlpha (initial)");
                isHeter = false;
                OptimizeHomoFixedAlpha(myMinimizer);
                PC[1] = PC[0];
                fn.globalPC2 = fn.globalPC;
                isHeter = true;
            }
            {
                PhaseTimer t("OptimizeHeterFixedAlpha");
                OptimizeHeterFixedAlpha(myMinimizer);
            }
        } else {
            std::cout << "Estimation from OptimizeHeter:" << std::endl;
            {
                PhaseTimer t("OptimizeHomo (initial)");
                isHeter = false;
                OptimizeHomo(myMinimizer);
                PC[1] = PC[0];
                fn.globalPC2 = fn.globalPC;
                isHeter = true;
            }
            {
                PhaseTimer t("OptimizeHeter");
                OptimizeHeter(myMinimizer);
            }
        }
        if (fn.globalAlpha >= 0.5) {
            std::swap(fn.globalPC[0], fn.globalPC2[0]);
            std::swap(fn.globalPC[1], fn.globalPC2[1]);
        }
    }

    {
        PhaseTimer t("Calculate null-model LLK");
        fn.CalculateLLK0();
    }
    std::cout << "Contaminating Sample ";
    for (int i = 0; i < numPC; ++i) {
        std::cout << "PC" << i + 1 << ":" << fn.globalPC[i] << "\t";
    }
    std::cout << std::endl;
    std::cout << "Intended Sample ";
    for (int i = 0; i < numPC; ++i) {
        std::cout << "PC" << i + 1 << ":" << fn.globalPC2[i] << "\t";
    }
    std::cout << std::endl;
    std::cout << "FREEMIX(Alpha):" << (fn.globalAlpha < 0.5 ? fn.globalAlpha : (1 - fn.globalAlpha)) << std::endl;

    std::string fileName(OutputPrefix + ".Ancestry");
    std::ofstream fout(fileName);
    if(not fout.is_open())
    {
      error("Open file %s failed!",fileName.c_str());
      exit(EXIT_FAILURE);
    }
    //std::cout << "PC\tContaminatingSample\tIntendedSample"<<std::endl;
    fout<< "PC\tContaminatingSample\tIntendedSample"<<std::endl;
    for (int i = 0; i < numPC; ++i) {
      //std::cout << i + 1 << "\t" << fn.globalPC[i] << "\t"<< fn.globalPC2[i] << std::endl;
      fout << i + 1 << "\t" << fn.globalPC[i] << "\t"<< fn.globalPC2[i] << std::endl;
    }

    fout.close();

    if(!fout)
    {
        error("Errors detected when writing to file %s !",fileName.c_str());
        exit(EXIT_FAILURE);
    }
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
    startingPoint[numPC * 2] = FullLLKFunc::Logit(alpha);

    if (verbose) {
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
    alpha = FullLLKFunc::InvLogit(myMinimizer.point[numPC * 2]);
    for (int i = 0; i < numPC; ++i) {
        PC[0][i] = myMinimizer.point[i];
    }
    for (int i = numPC; i < numPC * 2; ++i) {
        PC[1][i - numPC] = myMinimizer.point[i];
    }

    if (ret == std::numeric_limits<double>::max()) return false;
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

    if (verbose) {
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
    startingPoint[numPC] = FullLLKFunc::Logit(alpha);
    if (verbose) {
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
    alpha = FullLLKFunc::InvLogit(myMinimizer.point[numPC]);
    for (int i = 0; i < numPC; ++i) {
        PC[0][i] = myMinimizer.point[i];
    }
    if (ret == std::numeric_limits<double>::max()) return false;
    else return true;
}

bool ContaminationEstimator::OptimizeHomoFixedAlpha(AmoebaMinimizer &myMinimizer) {
    Vector startingPoint("TestPoint", numPC);
    for (int i = 0; i < numPC; ++i) {
        startingPoint[i] = PC[0][i];
    }
    if (verbose) {
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
    startingPoint[0] = FullLLKFunc::Logit(alpha);

    if (verbose) {
        std::cerr << "Start point";
        std::cerr << "alpha:\t" << alpha << std::endl;
    }

    startingPoint.label = "startPoint";
    myMinimizer.func = &fn;
    myMinimizer.Reset(1);
    myMinimizer.point = startingPoint;
    double ret = myMinimizer.Minimize(epsilon);
    alpha = FullLLKFunc::InvLogit(myMinimizer.point[0]);
    if (ret == std::numeric_limits<double>::max()) return false;//not converge
    else return true;
}

int ContaminationEstimator::ReadSVDMatrix(const std::string &UDpath, const std::string &PCpath, const std::string &Mean) {
    ReadMatrixUD(UDpath);
//    ReadMatrixPC(PCpath);
    ReadMean(Mean);

    return 0;
}

int ContaminationEstimator::ReadMatrixUD(const std::string &path) {
    InputFile fin;
    std::string line;
    std::vector<PCtype> tmpUD(numPC, 0);
    if (!fin.openFile(path.c_str(), "r", InputFile::ifileCompression::DEFAULT)) {
        std::cerr << "Open file:" << path << "\t failed, exit!";
        exit(EXIT_FAILURE);
    }
    while (fin.readLine(line)==0) {
        std::stringstream ss(line);
        int index = 0;
        while(index < numPC && ss>>tmpUD[index])
        {
          index++;
        }
        // Upon finish the line, index == numPC is expected
        if (index < numPC) {
          warning("--NumPC should be less than or equal to the number of PCs in SVD files provided by --SVDPrefix! (Expected:%d vs Observed:%d)", numPC, index);
          warning("--NumPC only permits as large as 4 PCs when using SVD files in ${verifybamID}/resource/ directory!");
          warning("You can always prepare you own SVD files with arbitrary number of PCs with --RefVCF enabled.");
          exit(EXIT_FAILURE);
        }
        UD.push_back(tmpUD);
        //Initialize arrays
        NumMarker++;
        AFs.push_back(0.);
      line="";
    }
    AF2s.assign(AFs.begin(), AFs.end());
    fin.ifclose();
    return 0;
}

int ContaminationEstimator::ReadMatrixPC(const std::string &path) {
    InputFile fin;
    std::string line;
    std::vector<PCtype> tmpPC(numPC, 0);
    std::vector<std::vector<PCtype> > PCvec;
    if (!fin.openFile(path.c_str(), "r", InputFile::ifileCompression::DEFAULT)) {
        std::cerr << "Open file:" << path << "\t failed, exit!";
        exit(EXIT_FAILURE);
    }
    std::string sampleID;
    while (fin.readLine(line)==0) {
        std::stringstream ss(line);
        ss >> sampleID;
        for (int index = 0; index != numPC; ++index)
            ss >> tmpPC[index];
        PCvec.push_back(tmpPC);
      line="";
    }
    fin.ifclose();

    // calculate the mean and variance of PCs
    std::vector<double> sumv(numPC,0.);
    std::vector<double> ssqv(numPC,0.);

    for(int32_t i=0; i < (int32_t)PCvec.size(); ++i) {//each instance
        for(int32_t j=0; j < numPC; ++j) {//each dimension
            sumv[j] += PCvec[i][j];
            ssqv[j] += PCvec[i][j] * PCvec[i][j];
        }
    }
    for(int32_t i=0; i < numPC; ++i) {
        muv[i] = sumv[i] / PCvec.size();
        sdv[i] = sqrt(ssqv[i]/PCvec.size() - muv[i] * muv[i]);
    }

    return 0;
}

int ContaminationEstimator::ReadChooseBed(const std::string &path) {
  InputFile fin;
  std::string line, chr;
    int index(0), pos(0);
    char ref(0), alt(0);

  if (!fin.openFile(path.c_str(), "r", InputFile::ifileCompression::DEFAULT)) {
        std::cerr << "Open file:" << path << "\t failed, exit!";
        exit(EXIT_FAILURE);
    }
  while (fin.readLine(line)==0) {
        index++;
        std::stringstream ss(line);
        //std::string chr;
        //int pos;
        ss >> chr >> pos >> pos;
        ss >> ref >> alt;

        BedVec.push_back(region_t(chr, pos - 1, pos));
        PosVec.push_back(make_pair(chr, pos));
        ChooseBed[chr][pos] = std::make_pair(ref, alt);
        line="";
  }
    fin.ifclose();
    return 0;
}

int ContaminationEstimator::ReadMean(const std::string &path) {
    InputFile fin;
    std::string line;
    double mu(0);
    std::string snpName, chr;
  if (!fin.openFile(path.c_str(), "r", InputFile::ifileCompression::DEFAULT)) {
        std::cerr << "Open file:" << path << "\t failed, exit!";
        exit(EXIT_FAILURE);
    }
  while (fin.readLine(line)==0) {
        std::stringstream ss(line);
        ss >> snpName;
        chr = snpName.substr(0, snpName.find(':', 0));
        ss >> mu;
        means.push_back(mu);
        line="";
  }
    fin.ifclose();
    return 0;
}

int ContaminationEstimator::ReadAF(const std::string &path) {
    std::ifstream fin(path);
    std::string line;
    uint32_t pos(0);
    double AF(0);
    std::string chr;
    char ref(0), alt(0);
//    int beg(0),end(0);
    if (!fin.is_open()) {
        std::cerr << "Open file:" << path << "\t failed, exit!";
        exit(EXIT_FAILURE);
    }
    while (std::getline(fin, line)) {
//        if(line[0]=='#'||line.find("INDEL")!=std::string::npos) continue;
        std::stringstream ss(line);
        ss >> chr;
        ss >> pos >> pos;
        ss >> ref >> alt;
//        ss>>snpName>>snpName>>snpName>>snpName>>snpName>>snpName;
//        beg=snpName.find("EUR_AF=");
//        beg+=7;
//        AF=atof(snpName.substr(beg,4).c_str());
        ss >> AF;
        knownAF[chr][pos] = AF;
    }
    return 0;
}

int ContaminationEstimator::ReadBam(const char *bamFile, const char *faiFile,
                                    const char *bedFile, mplp_conf_t *mplpPtr) {
  viewer = SimplePileupViewer(&BedVec, bamFile, faiFile, bedFile, mplpPtr, 1);
  return 0;
}

int ContaminationEstimator::ReadPileup(const std::string & pileupFile) {
    viewer = SimplePileupViewer(ChooseBed,pileupFile);
    isPileupInput = true;
    return 0;
}

//bool ContaminationEstimator::IsSanityCheckOK()
//{
//    int effectSite(0), tmpDepth(0);
//    std::vector<int> depthVec;
//    std::string chr;
//    int pos;
//    for (size_t i = 0; i < NumMarker; ++i) {
//        chr = PosVec[i].first;
//        pos = PosVec[i].second;
//        if (viewer.posIndex.find(chr) == viewer.posIndex.end()) {
//            continue;
//        } else if (viewer.posIndex[chr].find(pos) == viewer.posIndex[chr].end()) {
//            continue;
//        }
//        depthVec.push_back(viewer.GetBaseInfoAt(chr, pos).size());
//    }
//    std::sort(depthVec.begin(), depthVec.end());
//    viewer.firstQT = depthVec[NumMarker/4.0];
//    viewer.thirdQT = depthVec[NumMarker/4.0 * 3];
//
//    float range = viewer.thirdQT - viewer.firstQT;
//
//    viewer.firstQT -= 1.5 * range;
//    viewer.thirdQT += 1.5 * range;
//
//    effectSite=0;
//    for (size_t i = 0; i < NumMarker; ++i) {
//        chr = PosVec[i].first;
//        pos = PosVec[i].second;
//        tmpDepth = viewer.GetBaseInfoAt(chr, pos).size();
//        if(tmpDepth < viewer.firstQT||
//           tmpDepth > viewer.thirdQT ) continue;
//        effectSite++;
//    }
//
//    std::cerr<<"Mean Depth: "<<viewer.avgDepth<<"\n"
//             <<"First Quantile Depth: "<<viewer.firstQT<<"\n"
//             <<"Third Quantile Depth: "<<viewer.thirdQT<<std::endl;
//    std::cerr<<"Removed "<<NumMarker-effectSite<<" out of "<< NumMarker <<" outlier sites."<<std::endl;
//    return double(effectSite)/NumMarker > 0.5 and effectSite > 7000;
//}

bool ContaminationEstimator::IsSanityCheckOK()
{
//    if(isPileupInput) return true;
    int tmpDepth(0);
//    std::vector<int> depthVec;
    std::string chr;
    notice("Number of marker in Reference Matrix:%d", NumMarker);
    notice("Number of marker shared with input file:%d", viewer.GetNumMarker());

    int pos;
    for (size_t i = 0; i < NumMarker; ++i) {
        chr = PosVec[i].first;
        pos = PosVec[i].second;
        if (viewer.posIndex.find(chr) == viewer.posIndex.end()) {
            continue;
        } else if (viewer.posIndex[chr].find(pos) == viewer.posIndex[chr].end()) {
            continue;
        }
        tmpDepth = viewer.GetBaseInfoAt(chr, pos).size();
        viewer.sdDepth += tmpDepth * tmpDepth;
    }

    viewer.sdDepth = sqrt(viewer.sdDepth/viewer.effectiveNumSite - viewer.avgDepth * viewer.avgDepth);

    double filteredDepth(0), filteredsdDepth(0);
    viewer.effectiveNumSite=0;

    for (size_t i = 0; i < NumMarker; ++i) {
        chr = PosVec[i].first;
        pos = PosVec[i].second;
        if (viewer.posIndex.find(chr) == viewer.posIndex.end()) {
            continue;
        } else if (viewer.posIndex[chr].find(pos) == viewer.posIndex[chr].end()) {
            continue;
        }
        tmpDepth = viewer.GetBaseInfoAt(chr, pos).size();
        if(tmpDepth == 0 || tmpDepth < (viewer.avgDepth - 3 * viewer.sdDepth)||
           tmpDepth > (viewer.avgDepth + 3 * viewer.sdDepth)) continue;
        viewer.effectiveNumSite++;
    }
    notice("Mean Depth:%f",viewer.avgDepth);
    notice("SD Depth:%f",viewer.sdDepth);
    notice("%d SNP markers remained after sanity check.", viewer.GetNumMarker());
    return viewer.GetNumMarker() > 1000 and viewer.GetNumMarker() > (NumMarker * 0.1);
}