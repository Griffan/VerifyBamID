// Test helper (not shipped): turn an aligned BAM into a directional
// methyl-seq BAM by applying the bisulfite conversion that methylation mode is
// designed to handle. In genome-forward coordinates (the orientation htslib
// stores SEQ in), a read from the "top"/CT strand has every C read as T and a
// read from the "bottom"/GA strand has every G read as A. Converting *every*
// applicable base simulates fully unmethylated DNA -- the maximal-conversion
// case -- and is deterministic, so the resulting VerifyBamID --Methylation
// output can be diffed against a committed golden file.
//
// The strand rule below is written out explicitly here, independently of
// MethylationModel, so the end-to-end test genuinely exercises that code path
// rather than tautologically agreeing with it.
//
// Usage: MethylateTestBam <in.bam> <out.bam>   (writes out.bam and out.bam.bai)

#include "htslib/sam.h"

#include <cstdint>
#include <cstdio>

// Set the 4-bit base code at query position i in a packed SEQ array, matching
// htslib's bam_seqi() nibble layout (even i -> high nibble, odd i -> low).
static inline void set_base(uint8_t *seq, int i, int code) {
    int shift = (~i & 1) << 2;
    seq[i >> 1] = (seq[i >> 1] & ~(0xf << shift)) | (code << shift);
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: %s <in.bam> <out.bam>\n", argv[0]);
        return 2;
    }
    samFile *in = sam_open(argv[1], "r");
    if (!in) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    sam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr) { fprintf(stderr, "cannot read header of %s\n", argv[1]); return 1; }
    samFile *out = sam_open(argv[2], "wb");
    if (!out || sam_hdr_write(out, hdr) < 0) {
        fprintf(stderr, "cannot write %s\n", argv[2]); return 1;
    }

    // seq_nt16 base codes: A=1, C=2, G=4, T=8.
    const int A = 1, C = 2, G = 4, T = 8;
    bam1_t *b = bam_init1();
    int ret;
    while ((ret = sam_read1(in, hdr, b)) >= 0) {
        const uint16_t f = b->core.flag;
        const bool paired = f & BAM_FPAIRED;
        const bool rev    = f & BAM_FREVERSE;
        const bool read1  = f & BAM_FREAD1;
        const bool mrev   = f & BAM_FMREVERSE;
        // Directional convention: the conversion channel is the strand read 1
        // aligned to. Read 1 uses its own orientation; read 2 uses its mate's.
        const bool read1Reverse = paired ? (read1 ? rev : mrev) : rev;
        const bool top = !read1Reverse;   // top => CT-converted, else GA-converted

        uint8_t *seq = bam_get_seq(b);
        for (int i = 0; i < b->core.l_qseq; ++i) {
            int base = bam_seqi(seq, i);
            if (top && base == C)      set_base(seq, i, T);
            else if (!top && base == G) set_base(seq, i, A);
        }
        if (sam_write1(out, hdr, b) < 0) { fprintf(stderr, "write error\n"); return 1; }
    }
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(in);
    if (sam_close(out) < 0) { fprintf(stderr, "error closing %s\n", argv[2]); return 1; }
    if (ret < -1) { fprintf(stderr, "error reading %s\n", argv[1]); return 1; }
    if (sam_index_build(argv[2], 0) < 0) { fprintf(stderr, "cannot index %s\n", argv[2]); return 1; }
    return 0;
}
