/// Integration test for SVDcalculator::ReadVcf genotype parsing.
///
/// Reads a test VCF that encodes expected AF and N_MISSING in each row's INFO
/// field, calls ReadVcf, then validates the genotype matrix produces matching
/// values. This exercises PL/GL/GT parsing, per-sample fallback, partial-missing
/// handling, and phased genotypes.

#include "SVDcalculator.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Expected {
    std::string testDesc;
    double af;
    int nMissing;
};

/// Parse EXPECTED_AF, EXPECTED_N_MISSING, and TESTDESC from the INFO column of
/// each data line in the VCF. Only includes rows with FILTER=PASS (matching
/// ReadVcf's filtering behavior).
static std::vector<Expected> parseExpected(const std::string& vcfPath) {
    std::vector<Expected> result;
    std::ifstream fin(vcfPath);
    std::string line;

    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;

        // Split line by tab
        std::vector<std::string> cols;
        std::istringstream iss(line);
        std::string col;
        while (std::getline(iss, col, '\t')) cols.push_back(col);
        if (cols.size() < 8) continue;

        // Only include PASS rows (ReadVcf skips non-PASS)
        if (cols[6] != "PASS") continue;

        Expected ev;
        ev.testDesc = "unknown";
        ev.af = -1.0;
        ev.nMissing = -1;

        // Parse semicolon-delimited INFO field
        std::istringstream infoStream(cols[7]);
        std::string field;
        while (std::getline(infoStream, field, ';')) {
            size_t eq = field.find('=');
            if (eq == std::string::npos) continue;
            std::string key = field.substr(0, eq);
            std::string val = field.substr(eq + 1);
            if (key == "TESTDESC")           ev.testDesc = val;
            else if (key == "EXPECTED_AF")   ev.af = std::stod(val);
            else if (key == "EXPECTED_N_MISSING") ev.nMissing = std::stoi(val);
        }

        result.push_back(ev);
    }
    return result;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <test.vcf>" << std::endl;
        return 1;
    }
    std::string vcfPath = argv[1];

    // Run ReadVcf
    SVDcalculator calc;
    std::vector<std::vector<char>> genotype;
    int nSamples = 0, nMarkers = 0;
    std::unordered_set<std::string> includeChr;
    calc.ReadVcf(vcfPath, genotype, nSamples, nMarkers, includeChr);

    // Parse expected values from VCF INFO fields
    std::vector<Expected> expected = parseExpected(vcfPath);

    if (nMarkers != static_cast<int>(expected.size())) {
        std::cerr << "FAIL: expected " << expected.size()
                  << " markers but ReadVcf returned " << nMarkers << std::endl;
        return 1;
    }

    int failures = 0;
    for (int m = 0; m < nMarkers; m++) {
        const auto& geno = genotype[m];
        const auto& exp  = expected[m];

        // Compute observed N_MISSING and AF from the genotype vector
        int nMissing    = 0;
        int sumGeno     = 0;
        int nNonMissing = 0;
        for (int s = 0; s < nSamples; s++) {
            if (geno[s] < 0) {
                nMissing++;
            } else {
                sumGeno += geno[s];
                nNonMissing++;
            }
        }
        double af = (nNonMissing > 0)
            ? static_cast<double>(sumGeno) / (2.0 * nNonMissing)
            : 0.0;

        bool ok = true;
        if (nMissing != exp.nMissing) {
            std::cerr << "FAIL [" << exp.testDesc << "]: N_MISSING expected "
                      << exp.nMissing << " got " << nMissing << std::endl;
            ok = false;
        }
        if (std::fabs(af - exp.af) > 0.001) {
            std::cerr << "FAIL [" << exp.testDesc << "]: AF expected "
                      << exp.af << " got " << af << std::endl;
            ok = false;
        }

        if (ok) {
            std::cerr << "PASS [" << exp.testDesc << "]" << std::endl;
        } else {
            failures++;
        }
    }

    if (failures > 0) {
        std::cerr << failures << " of " << nMarkers << " test(s) FAILED" << std::endl;
        return 1;
    }

    std::cerr << "All " << nMarkers << " tests passed" << std::endl;
    return 0;
}
