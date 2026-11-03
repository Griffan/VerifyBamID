// Generate a small coordinate-sorted BAM for the command-line read-filter test.
// One overlapping proper pair and one unpaired read all cover position 100.
// The default overlap correction must count the pair as one observation.
#include "htslib/faidx.h"
#include "htslib/sam.h"

#include <cstdio>
#include <cstdint>
#include <fstream>
#include <string>

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s <output-prefix>\n", argv[0]);
        return 2;
    }
    const std::string prefix = argv[1];
    const std::string reference = prefix + ".fa";
    const std::string sam = prefix + ".sam";
    const std::string bam = prefix + ".bam";
    // Avoid ambiguous alignments/BAQ in a homopolymer reference.
    std::string referenceSequence(300, 'A');
    uint32_t state = 12345;
    for (char &base : referenceSequence) {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        base = "ACGT"[state % 4];
    }
    referenceSequence[99] = 'A';
    {
        std::ofstream out(reference);
        out << ">chr1\n" << referenceSequence << "\n";
        if (!out) return 1;
    }
    if (fai_build(reference.c_str()) < 0) return 1;
    {
        std::ofstream out(sam);
        out << "@HD\tVN:1.6\tSO:coordinate\n"
            << "@SQ\tSN:chr1\tLN:300\n"
            << "@RG\tID:rg1\tSM:sample\n";
        const std::string quality(100, 'I');
        out << "pair\t99\tchr1\t41\t60\t100M\t=\t51\t110\t"
            << referenceSequence.substr(40, 100) << '\t' << quality << "\tRG:Z:rg1\n";
        out << "single\t0\tchr1\t41\t60\t100M\t*\t0\t0\t"
            << referenceSequence.substr(40, 100) << '\t' << quality << "\tRG:Z:rg1\n";
        out << "pair\t147\tchr1\t51\t60\t100M\t=\t41\t-110\t"
            << referenceSequence.substr(50, 100) << '\t' << quality << "\tRG:Z:rg1\n";
        if (!out) return 1;
    }
    samFile *in = sam_open(sam.c_str(), "r");
    if (!in) return 1;
    sam_hdr_t *header = sam_hdr_read(in);
    samFile *out = sam_open(bam.c_str(), "wb");
    if (!header || !out || sam_hdr_write(out, header) < 0) return 1;
    bam1_t *record = bam_init1();
    int status;
    while ((status = sam_read1(in, header, record)) >= 0) {
        if (sam_write1(out, header, record) < 0) return 1;
    }
    bam_destroy1(record);
    sam_hdr_destroy(header);
    sam_close(in);
    if (sam_close(out) < 0 || status < -1) return 1;
    if (sam_index_build(bam.c_str(), 0) < 0) return 1;
    return 0;
}
