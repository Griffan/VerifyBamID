// Generates the long-read pileup test data: a two-contig reference and the same
// coordinate-sorted long reads as an indexed BAM and an indexed CRAM.
//
// Every marker position used by TestLongReadPileup.cmake holds reference base A:
//   chr1:2001, chr1:2002  adjacent markers
//   chr1:4001             deleted in `del_alt`
//   chr1:15001            covered by no read
//   chr2:3001
// Reads (all MAPQ 60):
//   del_alt  chr1:1001 + 3000M1D1000M: base C (alt) at 2001, quality Q20 at 2002,
//            deletion over 4001, otherwise reference bases at Q40
//   rev      chr1:1501 - 3000M: reference bases at Q30
//   chr2     chr2:1001 + 4000M: reference bases at Q40
#include "htslib/faidx.h"
#include "htslib/sam.h"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

static const int kChr1Length = 20000;
static const int kChr2Length = 8000;

/// Pseudo-random sequence (xorshift), so alignments are unambiguous.
static std::string RandomSequence(int length, uint32_t seed) {
    std::string sequence(length, 'A');
    for (char &base : sequence) {
        seed ^= seed << 13;
        seed ^= seed >> 17;
        seed ^= seed << 5;
        base = "ACGT"[seed % 4];
    }
    return sequence;
}

/// Reference bases at 1-based positions [start, end].
static std::string Bases(const std::string &contig, int start, int end) {
    return contig.substr(start - 1, end - start + 1);
}

static std::string Record(const std::string &name, int flag, const std::string &contig, int pos,
                          const std::string &cigar, const std::string &bases, const std::string &quals) {
    return name + "\t" + std::to_string(flag) + "\t" + contig + "\t" + std::to_string(pos) + "\t60\t" + cigar +
           "\t*\t0\t0\t" + bases + "\t" + quals + "\tRG:Z:rg1\n";
}

/// Copies `sam` to `out` in the given htslib write mode ("wb" or "wc") and indexes it.
static bool Convert(const std::string &sam, const std::string &out, const char *mode, const std::string &reference) {
    samFile *in = sam_open(sam.c_str(), "r");
    samFile *writer = sam_open(out.c_str(), mode);
    if (!in || !writer || hts_set_fai_filename(writer, reference.c_str()) != 0) return false;
    sam_hdr_t *header = sam_hdr_read(in);
    if (!header || sam_hdr_write(writer, header) < 0) return false;
    bam1_t *record = bam_init1();
    int status;
    while ((status = sam_read1(in, header, record)) >= 0) {
        if (sam_write1(writer, header, record) < 0) return false;
    }
    bam_destroy1(record);
    sam_hdr_destroy(header);
    sam_close(in);
    if (sam_close(writer) < 0 || status < -1) return false;
    return sam_index_build(out.c_str(), 0) == 0;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s <output-prefix>\n", argv[0]);
        return 2;
    }
    const std::string prefix = argv[1];
    const std::string reference = prefix + ".fa";
    const std::string sam = prefix + ".sam";

    std::string chr1 = RandomSequence(kChr1Length, 12345);
    std::string chr2 = RandomSequence(kChr2Length, 67890);
    for (int pos : {2001, 2002, 4001, 15001}) chr1[pos - 1] = 'A';
    chr2[3001 - 1] = 'A';
    {
        std::ofstream out(reference);
        out << ">chr1\n" << chr1 << "\n>chr2\n" << chr2 << "\n";
        if (!out) return 1;
    }
    if (fai_build(reference.c_str()) < 0) return 1;

    std::string delAltBases = Bases(chr1, 1001, 4000) + Bases(chr1, 4002, 5001);
    std::string delAltQuals(delAltBases.size(), 'I');
    delAltBases[2001 - 1001] = 'C';
    delAltQuals[2002 - 1001] = '5';

    std::vector<std::string> records;
    records.push_back(Record("del_alt", 0, "chr1", 1001, "3000M1D1000M", delAltBases, delAltQuals));
    records.push_back(Record("rev", 16, "chr1", 1501, "3000M", Bases(chr1, 1501, 4500), std::string(3000, '?')));
    records.push_back(Record("chr2", 0, "chr2", 1001, "4000M", Bases(chr2, 1001, 5000), std::string(4000, 'I')));
    {
        std::ofstream out(sam);
        out << "@HD\tVN:1.6\tSO:coordinate\n"
            << "@SQ\tSN:chr1\tLN:" << kChr1Length << "\n"
            << "@SQ\tSN:chr2\tLN:" << kChr2Length << "\n"
            << "@RG\tID:rg1\tSM:sample\n";
        for (const std::string &record : records) out << record;
        if (!out) return 1;
    }

    if (!Convert(sam, prefix + ".bam", "wb", reference)) return 1;
    if (!Convert(sam, prefix + ".cram", "wc", reference)) return 1;

    // The first half of the BAM with the whole file's index, so reading the
    // reads fails part-way.
    std::ifstream bam(prefix + ".bam", std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(bam)), std::istreambuf_iterator<char>());
    std::ofstream truncated(prefix + ".truncated.bam", std::ios::binary);
    truncated.write(bytes.data(), bytes.size() / 2);
    std::ifstream index(prefix + ".bam.bai", std::ios::binary);
    std::ofstream truncatedIndex(prefix + ".truncated.bam.bai", std::ios::binary);
    truncatedIndex << index.rdbuf();
    return truncated && truncatedIndex ? 0 : 1;
}
