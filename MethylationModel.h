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

// Whether an observation from a read in channel `cs` can still identify one of
// the marker's two alleles unambiguously. The channel's conversion rewrites one
// base (C->T on Ct reads, G->A on Ga reads), so an allele carrying that base can
// no longer be told apart from its converted form. A marker is therefore
// unusable on a channel whenever one of its alleles IS the channel's convertible
// base -- regardless of what the other allele is:
//   - Ct  unusable whenever C is an allele (a genomic C may read as C or T);
//   - Ga  unusable whenever G is an allele (a genomic G may read as G or A).
// e.g. A/C is unusable on Ct even though the C never becomes the A. Such
// observations are dropped rather than collapsed back onto their allele -- the
// "no-collapse" selection rule. `ref`/`alt` are single-character alleles
// (case-insensitive); the Unknown channel is always unusable. Net effect: A/T
// usable on both strands, other transversions and transitions on one, C/G on
// neither.
bool observationUsable(char ref, char alt, ConversionStrand cs);

#endif // METHYLATIONMODEL_H_
