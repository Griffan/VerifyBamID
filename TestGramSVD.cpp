/// Test that the Gram matrix SVD approach (A^T*A eigendecomposition) produces
/// mathematically equivalent results to Eigen's JacobiSVD.
///
/// Generates a deterministic synthetic genotype matrix (500 markers x 50
/// samples, values in {0,1,2}), mean-centers it (matching ProcessRefVCF), then
/// computes the decomposition both ways and compares:
///   1. UD columns  (U * diag(sigma) from Jacobi  vs  A * V_k from Gram)
///   2. V columns   (right singular vectors)
///   3. Singular values  (direct from Jacobi  vs  sqrt(eigenvalues) from Gram)
///
/// SVD singular vectors have sign ambiguity -- both v_i and -v_i are valid
/// decompositions -- so each column comparison determines the sign alignment
/// via dot product before computing the error.
///
/// The test uses single-precision float (matching the production code's
/// MatrixXf), so we allow up to 1e-2 relative error per column -- empirically
/// the two methods agree to ~1e-5, well within this bound.

#include "Eigen/Dense"
#include <cmath>
#include <cstdint>
#include <cstdlib>
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

/// Compare the top K principal components from JacobiSVD and Gram matrix
/// approaches on a given mean-centered matrix.  Returns the number of
/// failures (0 = all comparisons passed).
int compareJacobiVsGram(const char* label, const MatrixXf& A, int K,
                        float tolerance) {
    // Method 1: JacobiSVD (existing approach)
    JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
    MatrixXf jacobiUD = svd.matrixU() * svd.singularValues().asDiagonal();
    MatrixXf jacobiV  = svd.matrixV();

    // Method 2: Gram matrix approach
    MatrixXf G = A.transpose() * A;
    SelfAdjointEigenSolver<MatrixXf> eigSolver(G);
    VectorXf eigenvalues  = eigSolver.eigenvalues().reverse();
    MatrixXf eigenvectors = eigSolver.eigenvectors().rowwise().reverse();

    MatrixXf gramV_k  = eigenvectors.leftCols(K);
    MatrixXf gramUD_k = A * gramV_k;

    int failures = 0;

    // Compare UD columns (accounting for sign ambiguity)
    for (int pc = 0; pc < K; ++pc) {
        float dot = jacobiUD.col(pc).dot(gramUD_k.col(pc));
        float sign = (dot >= 0.0f) ? 1.0f : -1.0f;

        float colNorm = jacobiUD.col(pc).norm();
        float maxDiff = (jacobiUD.col(pc) - sign * gramUD_k.col(pc)).cwiseAbs().maxCoeff();
        float relErr  = (colNorm > 0.0f) ? maxDiff / colNorm : maxDiff;

        if (relErr > tolerance) {
            std::cerr << "FAIL [" << label << "]: UD column " << pc
                      << " relative error = " << relErr
                      << " (max abs diff = " << maxDiff << ")" << std::endl;
            failures++;
        } else {
            std::cerr << "PASS [" << label << "]: UD column " << pc
                      << " relative error = " << relErr << std::endl;
        }
    }

    // Compare V columns
    for (int pc = 0; pc < K; ++pc) {
        float dot = jacobiV.col(pc).dot(gramV_k.col(pc));
        float sign = (dot >= 0.0f) ? 1.0f : -1.0f;

        float colNorm = jacobiV.col(pc).norm();
        float maxDiff = (jacobiV.col(pc) - sign * gramV_k.col(pc)).cwiseAbs().maxCoeff();
        float relErr  = (colNorm > 0.0f) ? maxDiff / colNorm : maxDiff;

        if (relErr > tolerance) {
            std::cerr << "FAIL [" << label << "]: V column " << pc
                      << " relative error = " << relErr
                      << " (max abs diff = " << maxDiff << ")" << std::endl;
            failures++;
        } else {
            std::cerr << "PASS [" << label << "]: V column " << pc
                      << " relative error = " << relErr << std::endl;
        }
    }

    // Verify singular values match eigenvalues
    for (int pc = 0; pc < K; ++pc) {
        float jacobiSV = svd.singularValues()(pc);
        float gramSV   = std::sqrt(std::max(0.0f, eigenvalues(pc)));
        float relErr   = (jacobiSV > 0.0f)
            ? std::fabs(jacobiSV - gramSV) / jacobiSV
            : std::fabs(gramSV);

        if (relErr > tolerance) {
            std::cerr << "FAIL [" << label << "]: singular value " << pc
                      << " Jacobi=" << jacobiSV << " Gram=" << gramSV
                      << " relative error=" << relErr << std::endl;
            failures++;
        } else {
            std::cerr << "PASS [" << label << "]: singular value " << pc
                      << " Jacobi=" << jacobiSV << " Gram=" << gramSV << std::endl;
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

    if (totalFailures > 0) {
        std::cerr << totalFailures << " total comparison(s) FAILED" << std::endl;
        return 1;
    }

    std::cerr << "All Gram vs Jacobi comparisons passed across all test cases"
              << std::endl;
    return 0;
}
