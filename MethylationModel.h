#ifndef METHYLATIONMODEL_H_
#define METHYLATIONMODEL_H_

#include "htslib/sam.h"
#include <cstdint>

// Support for genotyping SNP markers from directional bisulfite (WGBS) and
// enzymatic methyl-seq (EM-seq) reads, where unmethylated C is read as T.
//
// In genome-forward coordinates (the orientation htslib stores) every read
// falls into one of two conversion channels:
//   - Ct: genomic C may read as C or T (the "top"/OT+CTOT channel);
//   - Ga: genomic G may read as G or A (the "bottom"/OB+CTOB channel).
// See docs/methylation-design.md for the full model and derivation.
enum class ConversionStrand : uint8_t { Ct, Ga, Unknown };

// Conversion channel of a read inferred from its FLAG alone, assuming a
// directional library. The channel is the strand read 1 aligned to:
//   - a paired read must be a proper pair (BAM_FPROPER_PAIR), so that
//     BAM_FMREVERSE reliably gives the mate's (i.e. R1's) orientation for R2;
//     any non-proper pair returns Unknown rather than guessing;
//   - an unpaired read uses its own alignment strand.
// Pure function of the flag bits so it can be tested exhaustively.
ConversionStrand strandFromFlag(uint16_t flag);

// Conversion channel of a read, preferring an aligner methylation tag over the
// FLAG inference. Resolution order (BISCUIT's get_bsstrand precedence):
//   1. YD  — bwa-meth (type Z) / BISCUIT (type A): f->Ct, r->Ga, u/c->Unknown
//   2. ZS  — BSMAP: leading '+' -> Ct, '-' -> Ga
//   3. XG  — Bismark/DRAGEN/BSBolt: "CT" -> Ct, "GA" -> Ga
//   4. strandFromFlag(b->core.flag)
// A tag that is present but carries an unrecognized value falls through to the
// next source; a YD explicitly marking the strand unknown ('u'/'c') returns
// Unknown without falling through, since the aligner has declared it unclear.
ConversionStrand conversionStrandOf(const bam1_t *b);

// Whether an observation from a read in channel `cs` unambiguously identifies
// one of the marker's two alleles. An observation is unusable exactly when the
// channel's conversion could turn one allele into the other:
//   - Ct  is ambiguous whenever C is an allele (a T could be a converted C);
//   - Ga  is ambiguous whenever G is an allele (an A could be a converted G).
// `ref`/`alt` are single-character alleles (case-insensitive). Unknown channel
// is always unusable. This is the no-collapse selection rule: markers keep only
// the strand(s) on which they stay unambiguous — A/T both, other transversions
// and transitions one, C/G neither.
bool observationUsable(char ref, char alt, ConversionStrand cs);

#endif // METHYLATIONMODEL_H_
