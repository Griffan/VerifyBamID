// Exercise the production likelihood with synthetic resolved markers, without
// requiring a BAM or an ancestry panel. The reference evaluates a closed-form
// count likelihood in long double, without the production per-read lookup table.
#include "ContaminationEstimator.h"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace {
int failures = 0;
int checks = 0;

void expectNear(double got, long double want, const char* label) {
    ++checks;
    if (!std::isfinite(got) || std::fabs(got - want) > 1e-9L * (1 + std::fabs(want))) {
        std::cerr << "FAIL " << label << ": got " << got << ", expected " << want << '\n';
        ++failures;
    }
}

long double reference(int refs, int alts, int others, int q, double af, double alpha) {
    af = std::max(0.00005, std::min(0.99995, af));
    const long double f = af;
    const long double priors[] = {(1-f)*(1-f), 2*f*(1-f), f*f};
    const long double error = std::pow(10.L, -q / 10.L);
    long double terms[9];
    long double largest = -std::numeric_limits<long double>::infinity();
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            const long double refFraction = alpha*(1-a/2.L) + (1-alpha)*(1-b/2.L);
            const long double pRef = refFraction*(1-error) + (1-refFraction)*error/3;
            const long double pAlt = (1-refFraction)*(1-error) + refFraction*error/3;
            long double term = std::log(priors[a]*priors[b]);
            if (refs) term += refs*std::log(pRef);
            if (alts) term += alts*std::log(pAlt);
            if (others) term += others*std::log(2*error/3);
            terms[a*3+b] = term;
            largest = std::max(largest, term);
        }
    }
    long double total = 0;
    for (long double term : terms) total += std::exp(term-largest);
    return largest + std::log(total);
}

void setMarker(ContaminationEstimator& c, int refs, int alts, int others, int q, double af) {
    std::vector<char> bases(refs, '.');
    bases.insert(bases.end(), alts, 'C');
    bases.insert(bases.end(), others, 'G');
    c.NumMarker = 1;
    c.numThread = 1;
    c.isAFknown = true;
    c.isSanityCheckDisabled = false;
    c.viewer.avgDepth = bases.size();
    c.viewer.sdDepth = 0;
    c.viewer.baseInfo = {bases};
    c.viewer.qualInfo = {std::vector<char>(bases.size(), static_cast<char>(q+33))};
    c.resolvedMarkers = {{0, 'C', af}};
    c.AFs.resize(1);
    c.AF2s.resize(1);
    c.fn.ptr = &c;
}
}

int main() {
    ContaminationEstimator c;
    // Ordinary and underflowing depths, zero/high quality, alpha endpoints,
    // and clamped allele frequencies. Each marker passes the depth sanity check.
    for (int depth : {2, 100, 1000, 1200, 8000}) {
        for (int q : {0, 10, 40, 93}) {
            for (double af : {0., 0.01, 0.5, 0.99, 1.}) {
                for (double alpha : {0., 0.01, 0.1, 0.5, 1.}) {
                    setMarker(c, depth/2, depth/2, 0, q, af);
                    expectNear(c.fn.ComputeMixLLKs({}, {}, alpha),
                               reference(depth/2, depth/2, 0, q, af, alpha), "balanced marker");
                }
            }
        }
    }
    for (int q : {0, 20, 93}) {
        setMarker(c, 0, 0, 8000, q, 0.5);
        expectNear(c.fn.ComputeMixLLKs({}, {}, 0.1),
                   8000*(std::log(2.L/3) - q/10.L*std::log(10.L)), "other bases");
    }
    setMarker(c, 13, 7, 3, 20, 0.2);
    expectNear(c.fn.ComputeMixLLKs({}, {}, 0.1), reference(13, 7, 3, 20, 0.2, 0.1), "mixed bases");
    setMarker(c, 0, 0, 0, 40, 0.5);
    expectNear(c.fn.ComputeMixLLKs({}, {}, 0.1), 0, "empty marker");
    setMarker(c, 600, 600, 0, 40, 0.5);
    c.resolvedMarkers[0].baseInfoIndex = -1;
    expectNear(c.fn.ComputeMixLLKs({}, {}, 0.1), 0, "absent marker");
    c.resolvedMarkers[0].baseInfoIndex = 0;
    c.viewer.avgDepth = 20;
    expectNear(c.fn.ComputeMixLLKs({}, {}, 0.1), 0, "depth outlier");
    c.isSanityCheckDisabled = true;
    expectNear(c.fn.ComputeMixLLKs({}, {}, 0.1), reference(600, 600, 0, 40, .5, .1), "sanity disabled");
    c.NumMarker = 2;
    c.resolvedMarkers.push_back(c.resolvedMarkers[0]);
    c.AFs.resize(2);
    c.AF2s.resize(2);
    c.numThread = 2;
    expectNear(c.fn.ComputeMixLLKs({}, {}, .1), 2*reference(600, 600, 0, 40, .5, .1), "marker sum");
    setMarker(c, 1, 0, 0, 0, 0);
    c.fn.min_af = c.fn.max_af = 0;
    const double impossible = c.fn.ComputeMixLLKs({}, {}, 0);
    ++checks;
    if (!std::isinf(impossible) || impossible >= 0) {
        std::cerr << "FAIL impossible marker must have negative-infinite likelihood\n";
        ++failures;
    }
    std::cout << checks << " likelihood checks, " << failures << " failures\n";
    return failures ? 1 : 0;
}
