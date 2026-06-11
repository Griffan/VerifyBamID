#ifndef VERIFYBAMID_SVDCALCULATOR_H
#define VERIFYBAMID_SVDCALCULATOR_H


#include <vector>
#include <map>
#include <unordered_set>
#include "SimplePileupViewer.h"
#include "Eigen/Dense"

class SVDcalculator {
private:
    int numIndividual;
    int numMarker;
#define PCtype double

    std::vector<std::vector<PCtype> > UD;
    std::vector<std::vector<PCtype> > PC;
    std::vector<std::string> Samples;
    std::vector<PCtype> Mu;
    BED chooseBed;
    std::vector<region_t> BedVec;
public:
    SVDcalculator();
    ~SVDcalculator();

    /// Read a reference VCF, compute SVD, and write .UD, .mu, .bed, .V files.
    /// @param VcfPath  path to the reference panel VCF
    /// @param includeChr  if non-empty, only markers on these chromosomes are used
    /// @param skipMinSampleCountCheck  if true, < 1000 samples becomes a warning not an error
    /// @param numSVDPCs  number of PCs to write (0 = all)
    /// @param useGramSVD  if true, use Gram matrix (A^T*A) eigendecomposition instead
    ///     of JacobiSVD.  This is mathematically equivalent but dramatically reduces
    ///     peak memory for tall-skinny matrices (M markers >> N samples) because only
    ///     an N*N Gram matrix and M*K result columns are allocated, rather than full
    ///     M*N internal SVD matrices.  Note this is only a win when M >> N; for panels
    ///     with large N the N*N Gram matrix and its eigendecomposition dominate.  The
    ///     A^T*A multiplication also benefits from OpenMP parallelization, which Eigen
    ///     applies automatically to the matrix product.
    void ProcessRefVCF(const std::string& VcfPath,
                       const std::unordered_set<std::string>& includeChr,
                       bool skipMinSampleCountCheck = false,
                       int numSVDPCs = 10,
                       bool useGramSVD = false);

    /// Parse a VCF into a genotype matrix.
    /// @param includeChr  if non-empty, only markers on these chromosomes are used
    int ReadVcf(const std::string &VcfPath,
                std::vector<std::vector<char> >& genotype,
                int & nSamples, int& nMarkers,
                const std::unordered_set<std::string>& includeChr);

    /// Decompose a mean-centered M x N genotype matrix into its top `numPCs`
    /// principal components via the Gram matrix (A^T*A) eigendecomposition.
    /// Mathematically equivalent to ComputeSvdJacobi but avoids materializing
    /// the M x N thin-U matrix; see ProcessRefVCF's @param useGramSVD for the
    /// memory/performance trade-offs.  Implemented as a free-standing static
    /// method so it can be unit-tested directly against ComputeSvdJacobi.
    /// @param centered      mean-centered genotype matrix (M markers x N samples)
    /// @param numPCs        number of principal components (columns) to return
    /// @param matrixUD[out] M x numPCs matrix, UD = A * V_k (== U * diag(sigma))
    /// @param matrixPC[out] N x numPCs matrix of right singular vectors (loadings)
    /// @param singularValues[out] all min(M,N) singular values, descending. The
    ///     full set (not just the top numPCs) is returned so callers can compute
    ///     the variance-explained denominator correctly.
    static void ComputeSvdGram(const Eigen::MatrixXf& centered, int numPCs,
                               Eigen::MatrixXf& matrixUD, Eigen::MatrixXf& matrixPC,
                               Eigen::VectorXf& singularValues);

    /// Same contract as ComputeSvdGram, but using Eigen's JacobiSVD directly.
    static void ComputeSvdJacobi(const Eigen::MatrixXf& centered, int numPCs,
                                 Eigen::MatrixXf& matrixUD, Eigen::MatrixXf& matrixPC,
                                 Eigen::VectorXf& singularValues);

    std::vector<std::vector<double>> GetUDMatrix();
    std::vector<std::vector<double>> GetPCMatrix();
    std::vector<PCtype> GetMuArray();
    BED GetchooseBed();
    std::vector<region_t> GetBedVec();
    /// @param numSVDPCs  number of PCs to write (0 = all available)
    void WriteSVD(const std::string &VcfPath, int numSVDPCs);
};


#endif //VERIFYBAMID_SVDCALCULATOR_H
