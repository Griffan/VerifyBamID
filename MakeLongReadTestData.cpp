// Generates the long-read pileup test data: a two-contig reference and the same
// coordinate-sorted long reads as an indexed BAM and an indexed CRAM.
//
// Every marker position used by TestLongReadPileup.cmake holds reference base A:
//   chr1:2001, chr1:2002    adjacent markers
//   chr1:4001               deleted in `del_alt`
//   chr1:4003               2 aligned bases after the deletion in `del_alt`
//   chr1:9998, chr1:9999    3 and 2 aligned bases before the deletion in `eqx`
//   chr1:12001, 13001,      each covered twice by one of the `fold_*` reads below,
//     14001, 17001          and once by a `plain_*` read
//   chr1:15001              covered by no read
//   chr1:21301, 24351       covered by `fold_tie`, and by `fold_three` and `plain_3`
//   chr1:27001              covered by two reads named `*`
//   chr1:28503, chr1:28504  2 and 3 aligned bases after the deletion in `eqx_left`
//   chr1:29598, chr1:29600  3 and 1 aligned bases before the insertion in `padded`
//   chr2:3001
//   chr2:5998, chr2:5999    3 and 2 aligned bases before the insertion in `near_ins`
//   chr2:6002, chr2:6003    2 and 3 aligned bases after it
//   chr2:7001, 7601, 7901   each covered twice by one of the `fold_*` reads below
//   chr2:9502               the first base after the reference skip in `skip_after_ins`
// Reads (all MAPQ 60). A fold-back read has a forward primary alignment and a
// reverse supplementary alignment over the same marker. Most primary
// alignments are at least 1000 bp, so the reads are detected as long.
//   del_alt             chr1:1001 + 3000M1D1000M: base C (alt) at 2001, quality Q20 at
//                       2002, deletion over 4001, otherwise reference bases at Q40
//   rev                 chr1:1501 - 3000M: reference bases at Q30
//   eqx                 chr1:9001 + 999=1X1D1000=: Q40; a mismatch at 10000, then a
//                       deletion over 10001
//   fold_near_indel     primary chr1:11501 500M1I500M at Q40, so its reference base at
//                       12001 follows the insertion; supplementary chr1:11901 200M at
//                       Q40 with base C (alt) at 12001
//   plain_1             chr1:11951 + 1100M: reference bases at Q20
//   fold_deletion       primary chr1:12501 500M1D500M at Q40, deleting 13001;
//                       supplementary chr1:12901 200M of reference bases at Q40
//   fold_best_filtered  primary chr1:13501 500M1I500M at Q40, so its reference base at
//                       14001 follows the insertion; supplementary chr1:13901 200M of
//                       reference bases at Q30
//   fold_low_quality    primary chr1:16501 1000M of reference bases at Q40;
//                       supplementary chr1:16951 100M at Q30 with base C (alt) at Q5
//                       at 17001
//   plain_2             chr1:16601 + 1000M: reference bases at Q20
//   fold_tie            primary chr1:21001 1000M; supplementary chr1:21201 200M; both
//                       reference bases at Q30
//   fold_three          primary chr1:24001 1000M and reverse supplementary chr1:24201
//                       200M of reference bases at Q40, and a forward supplementary
//                       chr1:24301 100M at Q30 with base C (alt) at 24351
//   plain_3             chr1:24101 + 1000M: reference bases at Q20
//   *, *                two reads without names, chr1:26501 + 1000M and chr1:26601 +
//                       1000M, the second with base C (alt) at 27001; all Q40
//   eqx_left            chr1:28001 + 500=1D1X500=: Q40; a deletion over 28501, then a
//                       mismatch at 28502
//   padded              chr1:29101 + 500M1P1I400M: reference bases at Q40
//   chr2                chr2:1001 + 4000M: reference bases at Q40
//   near_ins            chr2:5501 + 50S500M1I1000M10H: reference bases at Q40,
//                       insertion after 6000
//   fold_agree          primary chr2:6801 400M at Q30; supplementary chr2:6901 200M at
//                       Q40; both reference bases
//   fold_disagree       primary chr2:7401 400M at Q30 with base C (alt) at 7601;
//                       supplementary chr2:7501 200M of reference bases at Q40
//   plain               chr2:7551 + 200M: reference bases at Q20
//   fold_first          primary chr2:7851 1000M at Q40; supplementary chr2:7871 60M at
//                       Q30; both reference bases
//   skip_after_ins      chr2:8901 + 600M1I1N400M: reference bases at Q40, with an
//                       insertion directly before a 1-base reference skip
#include "htslib/faidx.h"
#include "htslib/sam.h"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

static const int kChr1Length = 32000;
static const int kChr2Length = 10000;

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
    for (int pos : {2001, 2002, 4001, 4003, 9998, 9999, 12001, 13001, 14001, 15001, 17001, 21301, 24351, 27001,
                    28503, 28504, 29598, 29600})
        chr1[pos - 1] = 'A';
    for (int pos : {3001, 5998, 5999, 6002, 6003, 7001, 7601, 7901, 9502}) chr2[pos - 1] = 'A';
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

    const int kSupplementaryReverse = 2048 | 16;
    std::vector<std::string> records;
    records.push_back(Record("del_alt", 0, "chr1", 1001, "3000M1D1000M", delAltBases, delAltQuals));
    records.push_back(Record("rev", 16, "chr1", 1501, "3000M", Bases(chr1, 1501, 4500), std::string(3000, '?')));

    const char mismatch = chr1[10000 - 1] == 'C' ? 'G' : 'C';
    records.push_back(Record("eqx", 0, "chr1", 9001, "999=1X1D1000=",
                             Bases(chr1, 9001, 9999) + mismatch + Bases(chr1, 10002, 11001), std::string(2000, 'I')));

    records.push_back(Record("fold_near_indel", 0, "chr1", 11501, "500M1I500M",
                             Bases(chr1, 11501, 12000) + "G" + Bases(chr1, 12001, 12500), std::string(1001, 'I')));
    std::string foldNearIndelSupplementary = Bases(chr1, 11901, 12100);
    foldNearIndelSupplementary[12001 - 11901] = 'C';
    records.push_back(Record("fold_near_indel", kSupplementaryReverse, "chr1", 11901, "200M",
                             foldNearIndelSupplementary, std::string(200, 'I')));
    records.push_back(Record("plain_1", 0, "chr1", 11951, "1100M", Bases(chr1, 11951, 13050), std::string(1100, '5')));

    records.push_back(Record("fold_deletion", 0, "chr1", 12501, "500M1D500M",
                             Bases(chr1, 12501, 13000) + Bases(chr1, 13002, 13501), std::string(1000, 'I')));
    records.push_back(Record("fold_deletion", kSupplementaryReverse, "chr1", 12901, "200M", Bases(chr1, 12901, 13100),
                             std::string(200, 'I')));

    records.push_back(Record("fold_best_filtered", 0, "chr1", 13501, "500M1I500M",
                             Bases(chr1, 13501, 14000) + "G" + Bases(chr1, 14001, 14500), std::string(1001, 'I')));
    records.push_back(Record("fold_best_filtered", kSupplementaryReverse, "chr1", 13901, "200M",
                             Bases(chr1, 13901, 14100), std::string(200, '?')));

    records.push_back(Record("fold_low_quality", 0, "chr1", 16501, "1000M", Bases(chr1, 16501, 17500),
                             std::string(1000, 'I')));
    records.push_back(Record("plain_2", 0, "chr1", 16601, "1000M", Bases(chr1, 16601, 17600), std::string(1000, '5')));
    std::string foldLowQualityBases = Bases(chr1, 16951, 17050);
    std::string foldLowQualityQuals(100, '?');
    foldLowQualityBases[17001 - 16951] = 'C';
    foldLowQualityQuals[17001 - 16951] = '&';
    records.push_back(Record("fold_low_quality", kSupplementaryReverse, "chr1", 16951, "100M", foldLowQualityBases,
                             foldLowQualityQuals));

    records.push_back(Record("fold_tie", 0, "chr1", 21001, "1000M", Bases(chr1, 21001, 22000), std::string(1000, '?')));
    records.push_back(Record("fold_tie", kSupplementaryReverse, "chr1", 21201, "200M", Bases(chr1, 21201, 21400),
                             std::string(200, '?')));

    records.push_back(Record("fold_three", 0, "chr1", 24001, "1000M", Bases(chr1, 24001, 25000), std::string(1000, 'I')));
    records.push_back(Record("plain_3", 0, "chr1", 24101, "1000M", Bases(chr1, 24101, 25100), std::string(1000, '5')));
    records.push_back(Record("fold_three", kSupplementaryReverse, "chr1", 24201, "200M", Bases(chr1, 24201, 24400),
                             std::string(200, 'I')));
    std::string foldThreeThirdBases = Bases(chr1, 24301, 24400);
    foldThreeThirdBases[24351 - 24301] = 'C';
    records.push_back(Record("fold_three", 2048, "chr1", 24301, "100M", foldThreeThirdBases, std::string(100, '?')));

    records.push_back(Record("*", 0, "chr1", 26501, "1000M", Bases(chr1, 26501, 27500), std::string(1000, 'I')));
    std::string unnamedAltBases = Bases(chr1, 26601, 27600);
    unnamedAltBases[27001 - 26601] = 'C';
    records.push_back(Record("*", 0, "chr1", 26601, "1000M", unnamedAltBases, std::string(1000, 'I')));

    const char leftMismatch = chr1[28502 - 1] == 'C' ? 'G' : 'C';
    records.push_back(Record("eqx_left", 0, "chr1", 28001, "500=1D1X500=",
                             Bases(chr1, 28001, 28500) + leftMismatch + Bases(chr1, 28503, 29002),
                             std::string(1001, 'I')));
    records.push_back(Record("padded", 0, "chr1", 29101, "500M1P1I400M",
                             Bases(chr1, 29101, 29600) + "G" + Bases(chr1, 29601, 30000), std::string(901, 'I')));

    records.push_back(Record("chr2", 0, "chr2", 1001, "4000M", Bases(chr2, 1001, 5000), std::string(4000, 'I')));
    const std::string nearInsBases = std::string(50, 'T') + Bases(chr2, 5501, 6000) + "G" + Bases(chr2, 6001, 7000);
    records.push_back(Record("near_ins", 0, "chr2", 5501, "50S500M1I1000M10H", nearInsBases,
                             std::string(nearInsBases.size(), 'I')));
    records.push_back(Record("fold_agree", 0, "chr2", 6801, "400M", Bases(chr2, 6801, 7200), std::string(400, '?')));
    records.push_back(Record("fold_agree", kSupplementaryReverse, "chr2", 6901, "200M", Bases(chr2, 6901, 7100),
                             std::string(200, 'I')));
    std::string foldDisagreeBases = Bases(chr2, 7401, 7800);
    foldDisagreeBases[7601 - 7401] = 'C';
    records.push_back(Record("fold_disagree", 0, "chr2", 7401, "400M", foldDisagreeBases, std::string(400, '?')));
    records.push_back(Record("fold_disagree", kSupplementaryReverse, "chr2", 7501, "200M", Bases(chr2, 7501, 7700),
                             std::string(200, 'I')));
    records.push_back(Record("plain", 0, "chr2", 7551, "200M", Bases(chr2, 7551, 7750), std::string(200, '5')));
    records.push_back(Record("fold_first", 0, "chr2", 7851, "1000M", Bases(chr2, 7851, 8850), std::string(1000, 'I')));
    records.push_back(Record("fold_first", kSupplementaryReverse, "chr2", 7871, "60M", Bases(chr2, 7871, 7930),
                             std::string(60, '?')));
    records.push_back(Record("skip_after_ins", 0, "chr2", 8901, "600M1I1N400M",
                             Bases(chr2, 8901, 9500) + "G" + Bases(chr2, 9502, 9901), std::string(1001, 'I')));
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
