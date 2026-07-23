#include "MethylationModel.h"

#include <cctype>
#include <cstring>

ConversionStrand strandFromFlag(uint16_t flag) {
    if (flag & BAM_FPAIRED) {
        // A non-proper pair gives no reliable read-1 orientation: for read 2 we
        // read the mate's strand from BAM_FMREVERSE, which is only meaningful
        // when the pair mapped properly. Refuse to guess.
        if (!(flag & BAM_FPROPER_PAIR)) return ConversionStrand::Unknown;
        // The conversion channel is the strand read 1 aligned to. For read 1
        // that is its own reverse bit; for read 2 it is the mate's reverse bit.
        bool read1Reverse = (flag & BAM_FREAD1) ? (flag & BAM_FREVERSE)
                                                : (flag & BAM_FMREVERSE);
        return read1Reverse ? ConversionStrand::Ga : ConversionStrand::Ct;
    }
    return (flag & BAM_FREVERSE) ? ConversionStrand::Ga : ConversionStrand::Ct;
}

ConversionStrand conversionStrandOf(const bam1_t *b) {
    // 1. YD — bwa-meth (type Z) and BISCUIT (type A). 'u'/'c' mean the aligner
    //    could not determine the strand, so we stop rather than fall through.
    if (uint8_t *yd = bam_aux_get(b, "YD")) {
        char v = 0;
        if (*yd == 'A') v = bam_aux2A(yd);
        else if (*yd == 'Z') { const char *s = bam_aux2Z(yd); if (s) v = s[0]; }
        if (v == 'f') return ConversionStrand::Ct;
        if (v == 'r') return ConversionStrand::Ga;
        if (v == 'u' || v == 'c') return ConversionStrand::Unknown;
        // Any other value is unexpected; fall through to the next source.
    }
    // 2. ZS — BSMAP; the strand is the first character.
    if (uint8_t *zs = bam_aux_get(b, "ZS")) {
        if (*zs == 'Z') {
            const char *s = bam_aux2Z(zs);
            if (s && s[0] == '+') return ConversionStrand::Ct;
            if (s && s[0] == '-') return ConversionStrand::Ga;
        }
    }
    // 3. XG — Bismark / DRAGEN / BSBolt genome-conversion tag.
    if (uint8_t *xg = bam_aux_get(b, "XG")) {
        if (*xg == 'Z') {
            const char *s = bam_aux2Z(xg);
            if (s && strcmp(s, "CT") == 0) return ConversionStrand::Ct;
            if (s && strcmp(s, "GA") == 0) return ConversionStrand::Ga;
        }
    }
    // 4. No usable tag: infer from the FLAG.
    return strandFromFlag(b->core.flag);
}

bool observationUsable(char ref, char alt, ConversionStrand cs) {
    char r = static_cast<char>(std::toupper(static_cast<unsigned char>(ref)));
    char a = static_cast<char>(std::toupper(static_cast<unsigned char>(alt)));
    switch (cs) {
        // A T could be a converted C, so any marker with a C allele is ambiguous.
        case ConversionStrand::Ct: return r != 'C' && a != 'C';
        // An A could be a converted G, so any marker with a G allele is ambiguous.
        case ConversionStrand::Ga: return r != 'G' && a != 'G';
        case ConversionStrand::Unknown: return false;
    }
    return false;
}
