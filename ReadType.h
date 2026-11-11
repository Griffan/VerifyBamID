#ifndef VERIFYBAMID_READTYPE_H
#define VERIFYBAMID_READTYPE_H

#include <cstdint>
#include <string>
#include <vector>

/// Kind of sequencing reads in the input, which selects pileup defaults suited
/// to them. `Auto` is only an input value; detection always resolves it to
/// `Short` or `Long`.
enum class ReadType { Auto, Short, Long };

/// Reads whose median length (bp) is at least this are treated as long reads.
/// Short-read platforms top out around 300 bp, while long-read libraries have
/// medians in the kilobases.
constexpr int32_t kLongReadMinMedianLength = 1000;

/// Number of reads sampled from the start of the file when detecting the read
/// type.
constexpr int kReadTypeSampleSize = 500;

/// Most records read while looking for the sample, so a file with few usable
/// reads near its start is not read to the end.
constexpr int kReadTypeMaxRecordsExamined = 20000;

/// Value of --adjust-MQ meaning "use the default for the read type": 40 for
/// short reads, 0 (disabled) for long reads.
constexpr int kAdjustMQByReadType = -1;

/// Value of --DiscardNearIndel meaning "use the default for the read type":
/// kLongReadDiscardNearIndel for long reads, 0 (off) for short reads.
constexpr int kDiscardNearIndelByReadType = -1;

/// Long-read default for --DiscardNearIndel. When an indel error in a long
/// read sits next to a SNP, the aligner often places it on the far side and
/// reads the SNP as the reference base, which looks like contamination. BAQ
/// handles this for short reads but is off for long reads.
constexpr int kLongReadDiscardNearIndel = 2;

/// The resolved read type and a human-readable reason for it, for logging.
struct ReadTypeDecision {
    ReadType type;
    std::string reason;
};

/// Parses "auto", "short" or "long", ignoring case. Returns false and leaves
/// `out` unchanged if `text` is none of these.
bool ParseReadType(const std::string &text, ReadType *out);

/// Returns "auto", "short" or "long".
const char *ReadTypeName(ReadType type);

/// Classifies reads as short or long.
///
/// `platforms` holds the PL value of every @RG header line ("" where PL is
/// absent); `readLengths` holds the lengths of the sampled reads. The reads are
/// long if every @RG declares a long-read platform (ONT or PACBIO); otherwise
/// they are long if the median read length is at least
/// kLongReadMinMedianLength. With no @RG platforms and no sampled reads, the
/// result is short.
ReadTypeDecision ClassifyReadType(const std::vector<std::string> &platforms,
                                  std::vector<int32_t> readLengths);

/// Detects the read type of a SAM/BAM/CRAM file from its @RG PL values and the
/// query lengths (from CIGAR) of up to kReadTypeSampleSize primary, mapped,
/// non-duplicate, QC-passing reads from the start of the file. `refPath` is the
/// reference FASTA for CRAM, or null. Calls error() if the file cannot be read.
ReadTypeDecision DetectReadType(const char *path, const char *refPath);

/// Whether BAQ should be used for a resolved read type (Short or Long). It is
/// off for long reads because its cost grows with read length times the
/// read's net indel imbalance, which is prohibitive for long reads.
bool ReadTypeUsesBaq(ReadType type);

/// Resolves the --adjust-MQ value for a resolved read type: an explicit
/// `requested` value is returned unchanged, and kAdjustMQByReadType becomes 40
/// for short reads and 0 for long reads.
int ResolveAdjustMQ(ReadType type, int requested);

/// Resolves the --DiscardNearIndel value for a resolved read type: an explicit
/// `requested` value is returned unchanged, and kDiscardNearIndelByReadType
/// becomes kLongReadDiscardNearIndel for long reads and 0 for short reads.
int ResolveDiscardNearIndel(ReadType type, int requested);

/// Whether a resolved read type counts each read at most once per marker. Long
/// reads can align twice over one position under one name, most often an ONT
/// fold-back read whose supplementary alignment covers its own primary, and
/// both copies observe the same molecule. Short reads keep htslib's handling
/// of overlapping mates.
bool ReadTypeCountsEachReadOnce(ReadType type);

#endif // VERIFYBAMID_READTYPE_H
