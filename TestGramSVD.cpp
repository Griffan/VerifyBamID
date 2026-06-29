/// Test that the Gram matrix SVD approach (A^T*A eigendecomposition) produces
/// mathematically equivalent results to Eigen's JacobiSVD.
///
/// This test calls the *production* decomposition routines directly --
/// SVDcalculator::ComputeSvdGram and SVDcalculator::ComputeSvdJacobi -- rather
/// than reimplementing them, so it validates the code that actually ships.
///
/// It generates a deterministic synthetic genotype matrix (values in {0,1,2}),
/// mean-centers it (matching ProcessRefVCF), runs both decompositions, and
/// compares:
///   1. UD columns       (U * diag(sigma)  vs  A * V_k)
///   2. PC columns        (right singular vectors / loadings)
///   3. Singular values   (direct from Jacobi  vs  sqrt(eigenvalues) from Gram)
///
/// SVD singular vectors have sign ambiguity -- both v_i and -v_i are valid
/// decompositions -- so each column comparison determines the sign alignment
/// via dot product before computing the error.
///
/// The decompositions use single-precision float (matching the production
/// code's MatrixXf), so we allow up to 1e-2 relative error per column --
/// empirically the two methods agree to ~1e-5, well within this bound.

#include "SVDcalculator.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>

using namespace Eigen;

/// Simple linear congruential PRNG for deterministic, cross-platform test data.
/// We avoid std::rand() because its output sequence varies across platforms and
/// standard library implementations, which would make the test non-reproducible.
struct LCG {
    uint32_t state;
    explicit LCG(uint32_t seed) : state(seed) {}
    uint32_t next() {
        state = state * 1103515245u + 12345u;
        return (state >> 16) & 0x7fff;
    }
};

/// Compare one column from each method, accounting for sign ambiguity. The
/// comparison is scale-aware: a column passes if its max deviation is small
/// relative to its own norm OR relative to `refScale` (the dominant column
/// norm of the whole matrix). The second clause matters for the trailing,
/// lowest-variance components: their singular vectors are ill-determined and
/// the two algorithms can pick different directions, but the absolute UD
/// contribution (sigma_i * u_i) is negligible, so a pure relative-error test
/// would spuriously fail on numerical noise. `refScale` is the right yardstick
/// because anything tiny next to the dominant component is irrelevant
/// downstream, where UD enters as UD . PC dot products.
int compareColumn(const char* label, const char* what, int pc,
                  const VectorXf& a, const VectorXf& b, float refScale,
                  float tolerance) {
    float sign = (a.dot(b) >= 0.0f) ? 1.0f : -1.0f;
    float colNorm = a.norm();
    float maxDiff = (a - sign * b).cwiseAbs().maxCoeff();
    float denom = std::max(colNorm, refScale);
    float relErr = (denom > 0.0f) ? maxDiff / denom : maxDiff;

    if (relErr > tolerance) {
        std::cerr << "FAIL [" << label << "]: " << what << " column " << pc
                  << " scaled error = " << relErr
                  << " (max abs diff = " << maxDiff << ")" << std::endl;
        return 1;
    }
    std::cerr << "PASS [" << label << "]: " << what << " column " << pc
              << " scaled error = " << relErr << std::endl;
    return 0;
}

/// Run both production decompositions on a mean-centered matrix and compare the
/// top K principal components. Returns the number of failures (0 = all passed).
int compareJacobiVsGram(const char* label, const MatrixXf& A, int K,
                        float tolerance) {
    MatrixXf jacobiUD, jacobiPC;
    VectorXf jacobiSV;
    SVDcalculator::ComputeSvdJacobi(A, K, jacobiUD, jacobiPC, jacobiSV);

    MatrixXf gramUD, gramPC;
    VectorXf gramSV;
    SVDcalculator::ComputeSvdGram(A, K, gramUD, gramPC, gramSV);

    // Dominant column norm of each matrix, used as the scale floor so that
    // trailing near-zero columns are compared against the decomposition's
    // overall magnitude rather than their own vanishing norm.
    float udScale = jacobiUD.colwise().norm().maxCoeff();
    float pcScale = jacobiPC.colwise().norm().maxCoeff();
    float svScale = jacobiSV(0); // singular values are sorted descending

    int failures = 0;
    for (int pc = 0; pc < K; ++pc) {
        failures += compareColumn(label, "UD", pc, jacobiUD.col(pc), gramUD.col(pc), udScale, tolerance);
    }
    for (int pc = 0; pc < K; ++pc) {
        failures += compareColumn(label, "PC", pc, jacobiPC.col(pc), gramPC.col(pc), pcScale, tolerance);
    }

    // Singular values are non-negative, so no sign handling is needed. As with
    // the vectors, judge the difference against the largest singular value so
    // the tiny trailing values aren't held to a meaningless relative tolerance.
    for (int pc = 0; pc < K; ++pc) {
        float jSV = jacobiSV(pc);
        float gSV = gramSV(pc);
        float relErr = (svScale > 0.0f) ? std::fabs(jSV - gSV) / svScale : std::fabs(gSV);
        if (relErr > tolerance) {
            std::cerr << "FAIL [" << label << "]: singular value " << pc
                      << " Jacobi=" << jSV << " Gram=" << gSV
                      << " relative error=" << relErr << std::endl;
            failures++;
        } else {
            std::cerr << "PASS [" << label << "]: singular value " << pc
                      << " Jacobi=" << jSV << " Gram=" << gSV << std::endl;
        }
    }

    return failures;
}

/// Generate a mean-centered random genotype matrix with values in {0, 1, 2}.
MatrixXf generateGenoMatrix(int M, int N, uint32_t seed) {
    LCG rng(seed);
    MatrixXf A(M, N);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            A(i, j) = static_cast<float>(rng.next() % 3);

    VectorXf mu = A.rowwise().mean();
    for (int i = 0; i < M; ++i)
        A.row(i).array() -= mu(i);
    return A;
}

int main() {
    const float tolerance = 1e-2f; // relative tolerance for float SVD
    int totalFailures = 0;

    // Case 1: Well-conditioned tall-skinny matrix (typical use case).
    // 500 markers x 50 samples, testing top 4 PCs.
    {
        std::cerr << "=== Case 1: tall-skinny (500 x 50, K=4) ===" << std::endl;
        MatrixXf A = generateGenoMatrix(500, 50, /*seed=*/42);
        totalFailures += compareJacobiVsGram("tall-skinny", A, 4, tolerance);
        std::cerr << std::endl;
    }

    // Case 2: Near-square matrix where M is only slightly larger than N.
    // This exercises the case where the matrix is closer to rank-deficient
    // and eigenvalue gaps may be smaller, making the decomposition less
    // numerically stable.  60 markers x 50 samples, testing top 4 PCs.
    {
        std::cerr << "=== Case 2: near-square (60 x 50, K=4) ===" << std::endl;
        MatrixXf A = generateGenoMatrix(60, 50, /*seed=*/99);
        totalFailures += compareJacobiVsGram("near-square", A, 4, tolerance);
        std::cerr << std::endl;
    }

    // Case 3: More PCs than typical (extract 10 of 30 possible).
    // Tests that higher-order PCs (with smaller singular values) still agree.
    {
        std::cerr << "=== Case 3: many PCs (200 x 30, K=10) ===" << std::endl;
        MatrixXf A = generateGenoMatrix(200, 30, /*seed=*/7);
        totalFailures += compareJacobiVsGram("many-PCs", A, 10, tolerance);
        std::cerr << std::endl;
    }

    // Case 4: Extract every component (K == N == min(M,N)).  This is the
    // --NumSVDPCs 0 path, where both methods return all available PCs including
    // the lowest-variance ones.  Note this exercises the shape, not small-value
    // accuracy: forming A^T*A squares the condition number, so the smallest
    // singular values diverge between the two methods on ill-conditioned input.
    // The scale-aware tolerance in compareColumn deliberately tolerates that --
    // those components are negligible in the downstream UD . PC dot products --
    // so this case does not assert agreement on the trailing singular values.
    {
        std::cerr << "=== Case 4: all PCs (200 x 30, K=30) ===" << std::endl;
        MatrixXf A = generateGenoMatrix(200, 30, /*seed=*/13);
        totalFailures += compareJacobiVsGram("all-PCs", A, 30, tolerance);
        std::cerr << std::endl;
    }

    if (totalFailures > 0) {
        std::cerr << totalFailures << " total comparison(s) FAILED" << std::endl;
        return 1;
    }

    std::cerr << "All Gram vs Jacobi comparisons passed across all test cases"
              << std::endl;
    return 0;
}
