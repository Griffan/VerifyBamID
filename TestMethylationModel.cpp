/// Unit tests for the methylation conversion model (MethylationModel.h).
///
/// Covers the three public functions exhaustively over their (tiny) domains:
///   1. strandFromFlag  — every relevant FLAG combination (pairing, proper-pair,
///      read 1/2, forward/reverse), which is the path used when no aligner
///      methylation tag is present (e.g. minibwa / bwa-meth output).
///   2. observationUsable — all six unordered SNP allele classes on each
///      conversion channel, against the truth table in docs/methylation-design.md,
///      plus allele-order and case symmetry.
///   3. conversionStrandOf — tag parsing for YD (types Z and A), ZS and XG, the
///      YD-unknown case, and that a tag overrides the FLAG inference.
///
/// Follows the framework-free pattern of TestGramSVD.cpp: accumulate failures,
/// print PASS/FAIL per check, return non-zero if any failed.

#include "MethylationModel.h"
#include "htslib/sam.h"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>

namespace {

// BAM FLAG bits, named for readability in the test cases below.
constexpr uint16_t PAIRED = 0x1, PROPER = 0x2, REVERSE = 0x10,
                   MREVERSE = 0x20, READ1 = 0x40, READ2 = 0x80;

const char *name(ConversionStrand cs) {
    switch (cs) {
        case ConversionStrand::Ct: return "Ct";
        case ConversionStrand::Ga: return "Ga";
        case ConversionStrand::Unknown: return "Unknown";
    }
    return "?";
}

int failures = 0;

void expectStrand(const char *what, ConversionStrand got, ConversionStrand want) {
    if (got != want) {
        std::cerr << "FAIL [" << what << "]: got " << name(got)
                  << ", want " << name(want) << std::endl;
        ++failures;
    } else {
        std::cerr << "PASS [" << what << "]: " << name(got) << std::endl;
    }
}

void expectUsable(const char *what, bool got, bool want) {
    if (got != want) {
        std::cerr << "FAIL [" << what << "]: got " << (got ? "usable" : "drop")
                  << ", want " << (want ? "usable" : "drop") << std::endl;
        ++failures;
    } else {
        std::cerr << "PASS [" << what << "]: " << (got ? "usable" : "drop") << std::endl;
    }
}

// Build a minimal bam1_t carrying just a FLAG and (optionally) one aux tag. A
// fresh record has zero-length name/cigar/seq, so an appended tag sits at the
// start of the aux region and bam_aux_get resolves it correctly.
bam1_t *makeRead(uint16_t flag) {
    bam1_t *b = bam_init1();
    b->core.flag = flag;
    return b;
}
void addZ(bam1_t *b, const char *tag, const char *val) {
    bam_aux_append(b, tag, 'Z', strlen(val) + 1, (const uint8_t *) val);
}
void addA(bam1_t *b, const char *tag, char val) {
    bam_aux_append(b, tag, 'A', 1, (const uint8_t *) &val);
}

void testStrandFromFlag() {
    std::cerr << "=== strandFromFlag ===" << std::endl;
    // Unpaired: own strand.
    expectStrand("SE forward", strandFromFlag(0), ConversionStrand::Ct);
    expectStrand("SE reverse", strandFromFlag(REVERSE), ConversionStrand::Ga);
    // Proper pair, read 1: own strand.
    expectStrand("proper R1 fwd", strandFromFlag(PAIRED | PROPER | READ1), ConversionStrand::Ct);
    expectStrand("proper R1 rev", strandFromFlag(PAIRED | PROPER | READ1 | REVERSE), ConversionStrand::Ga);
    // Proper pair, read 2: the MATE's strand (BAM_FMREVERSE), not its own.
    expectStrand("proper R2, mate fwd", strandFromFlag(PAIRED | PROPER | READ2), ConversionStrand::Ct);
    expectStrand("proper R2, mate rev", strandFromFlag(PAIRED | PROPER | READ2 | MREVERSE), ConversionStrand::Ga);
    // Read 2 forward-mapped but mate forward too -> still Ct (own REVERSE ignored).
    expectStrand("proper R2 own-rev, mate fwd",
                 strandFromFlag(PAIRED | PROPER | READ2 | REVERSE), ConversionStrand::Ct);
    // Non-proper pairs are always Unknown, regardless of read/strand bits.
    expectStrand("paired-not-proper R1", strandFromFlag(PAIRED | READ1), ConversionStrand::Unknown);
    expectStrand("paired-not-proper R2", strandFromFlag(PAIRED | READ2 | MREVERSE), ConversionStrand::Unknown);
    expectStrand("paired-not-proper rev", strandFromFlag(PAIRED | READ1 | REVERSE), ConversionStrand::Unknown);
}

void testObservationUsable() {
    std::cerr << "=== observationUsable ===" << std::endl;
    // (ref, alt, usable-on-Ct, usable-on-Ga) per the design-doc truth table.
    struct Case { char ref, alt; bool ct, ga; const char *label; };
    const Case cases[] = {
        {'A', 'T', true,  true,  "A/T"},
        {'A', 'G', true,  false, "A/G"},
        {'C', 'T', false, true,  "C/T"},
        {'A', 'C', false, true,  "A/C"},
        {'G', 'T', true,  false, "G/T"},
        {'C', 'G', false, false, "C/G"},
    };
    for (const Case &c : cases) {
        std::string ct = std::string(c.label) + " Ct";
        std::string ga = std::string(c.label) + " Ga";
        std::string un = std::string(c.label) + " Unknown";
        expectUsable(ct.c_str(), observationUsable(c.ref, c.alt, ConversionStrand::Ct), c.ct);
        expectUsable(ga.c_str(), observationUsable(c.ref, c.alt, ConversionStrand::Ga), c.ga);
        expectUsable(un.c_str(), observationUsable(c.ref, c.alt, ConversionStrand::Unknown), false);
    }
    // Allele order does not matter (unordered {ref,alt}).
    expectUsable("T/A Ct (== A/T)", observationUsable('T', 'A', ConversionStrand::Ct), true);
    expectUsable("G/A Ga (== A/G)", observationUsable('G', 'A', ConversionStrand::Ga), false);
    // Case-insensitive.
    expectUsable("c/t Ct lowercase", observationUsable('c', 't', ConversionStrand::Ct), false);
    expectUsable("a/g Ga lowercase", observationUsable('a', 'g', ConversionStrand::Ga), false);
}

void testConversionStrandOf() {
    std::cerr << "=== conversionStrandOf (tags) ===" << std::endl;
    // A flag that would infer Ga on its own, so a tag result proves the tag won.
    const uint16_t flagGa = REVERSE; // SE reverse -> Ga via FLAG

    {   bam1_t *b = makeRead(flagGa); addZ(b, "YD", "f");
        expectStrand("YD:Z:f beats flag", conversionStrandOf(b), ConversionStrand::Ct); bam_destroy1(b); }
    {   bam1_t *b = makeRead(0); addZ(b, "YD", "r");
        expectStrand("YD:Z:r", conversionStrandOf(b), ConversionStrand::Ga); bam_destroy1(b); }
    {   bam1_t *b = makeRead(flagGa); addA(b, "YD", 'f');
        expectStrand("YD:A:f (BISCUIT)", conversionStrandOf(b), ConversionStrand::Ct); bam_destroy1(b); }
    {   bam1_t *b = makeRead(0); addA(b, "YD", 'r');
        expectStrand("YD:A:r (BISCUIT)", conversionStrandOf(b), ConversionStrand::Ga); bam_destroy1(b); }
    {   bam1_t *b = makeRead(0); addA(b, "YD", 'u');
        expectStrand("YD:A:u -> Unknown", conversionStrandOf(b), ConversionStrand::Unknown); bam_destroy1(b); }
    {   bam1_t *b = makeRead(flagGa); addZ(b, "ZS", "+-");
        expectStrand("ZS:Z:+ beats flag", conversionStrandOf(b), ConversionStrand::Ct); bam_destroy1(b); }
    {   bam1_t *b = makeRead(0); addZ(b, "ZS", "--");
        expectStrand("ZS:Z:-", conversionStrandOf(b), ConversionStrand::Ga); bam_destroy1(b); }
    {   bam1_t *b = makeRead(flagGa); addZ(b, "XG", "CT");
        expectStrand("XG:Z:CT beats flag", conversionStrandOf(b), ConversionStrand::Ct); bam_destroy1(b); }
    {   bam1_t *b = makeRead(0); addZ(b, "XG", "GA");
        expectStrand("XG:Z:GA", conversionStrandOf(b), ConversionStrand::Ga); bam_destroy1(b); }
    // No tag -> FLAG inference (this is the minibwa / bwa-meth path).
    {   bam1_t *b = makeRead(PAIRED | PROPER | READ1);
        expectStrand("no tag -> FLAG (proper R1 fwd)", conversionStrandOf(b), ConversionStrand::Ct); bam_destroy1(b); }
}

} // namespace

int main() {
    testStrandFromFlag();
    testObservationUsable();
    testConversionStrandOf();
    if (failures > 0) {
        std::cerr << failures << " check(s) FAILED" << std::endl;
        return 1;
    }
    std::cerr << "All methylation-model checks passed" << std::endl;
    return 0;
}
