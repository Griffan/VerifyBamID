#include "ReadType.h"

#include <algorithm>
#include <cctype>
#include <cstdio>

#include "Error.h"
#include "htslib/sam.h"

static std::string ToUpper(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return text;
}

static bool IsLongReadPlatform(const std::string &platform) {
    const std::string upper = ToUpper(platform);
    return upper == "ONT" || upper == "PACBIO";
}

bool ParseReadType(const std::string &text, ReadType *out) {
    const std::string upper = ToUpper(text);
    if (upper == "AUTO") *out = ReadType::Auto;
    else if (upper == "SHORT") *out = ReadType::Short;
    else if (upper == "LONG") *out = ReadType::Long;
    else return false;
    return true;
}

const char *ReadTypeName(ReadType type) {
    switch (type) {
    case ReadType::Auto: return "auto";
    case ReadType::Short: return "short";
    case ReadType::Long: return "long";
    }
    return "unknown";
}

ReadTypeDecision ClassifyReadType(const std::vector<std::string> &platforms,
                                  std::vector<int32_t> readLengths) {
    if (!platforms.empty() &&
        std::all_of(platforms.begin(), platforms.end(), IsLongReadPlatform)) {
        std::vector<std::string> distinct;
        for (const std::string &platform : platforms) {
            const std::string upper = ToUpper(platform);
            if (std::find(distinct.begin(), distinct.end(), upper) == distinct.end()) distinct.push_back(upper);
        }
        std::string joined;
        for (const std::string &platform : distinct) joined += (joined.empty() ? "" : "/") + platform;
        return {ReadType::Long, "@RG PL is " + joined};
    }

    if (readLengths.empty())
        return {ReadType::Short, "no primary mapped reads to sample; assuming short reads"};

    const size_t mid = readLengths.size() / 2;
    std::nth_element(readLengths.begin(), readLengths.begin() + mid, readLengths.end());
    const int32_t median = readLengths[mid];
    const ReadType type = median >= kLongReadMinMedianLength ? ReadType::Long : ReadType::Short;

    char reason[160];
    snprintf(reason, sizeof(reason), "median length %d bp over the first %zu primary mapped reads (long if >= %d bp)",
             median, readLengths.size(), kLongReadMinMedianLength);
    return {type, reason};
}

ReadTypeDecision DetectReadType(const char *path, const char *refPath) {
    samFile *fp = sam_open(path, "r");
    if (!fp) error("Could not open %s to detect the read type", path);
    if (refPath && hts_set_fai_filename(fp, refPath) != 0)
        error("Could not use reference %s for %s", refPath, path);
    // Lengths come from the CIGAR, so CRAM needs neither sequences nor the
    // reference, and reads stored with SEQ '*' still count.
    if (hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_POS | SAM_CIGAR) != 0 ||
        hts_set_opt(fp, CRAM_OPT_DECODE_MD, 0) != 0)
        error("Could not set decoding options for %s", path);

    sam_hdr_t *header = sam_hdr_read(fp);
    if (!header) error("Could not read the header of %s", path);

    std::vector<std::string> platforms;
    kstring_t platform = KS_INITIALIZE;
    const int numReadGroups = sam_hdr_count_lines(header, "RG");
    for (int i = 0; i < numReadGroups; ++i) {
        platforms.push_back(sam_hdr_find_tag_pos(header, "RG", i, "PL", &platform) == 0 ? platform.s : "");
    }
    ks_free(&platform);

    std::vector<int32_t> readLengths;
    bam1_t *record = bam_init1();
    const uint16_t skipFlags = BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP;
    int examined = 0, ret = 0;
    while (readLengths.size() < static_cast<size_t>(kReadTypeSampleSize) &&
           examined++ < kReadTypeMaxRecordsExamined && (ret = sam_read1(fp, header, record)) >= 0) {
        if (record->core.flag & skipFlags) continue;
        const int64_t length = bam_cigar2qlen(record->core.n_cigar, bam_get_cigar(record));
        if (length > 0) readLengths.push_back(static_cast<int32_t>(std::min<int64_t>(length, INT32_MAX)));
    }
    bam_destroy1(record);
    if (ret < -1) error("Failed reading %s while detecting the read type", path);
    sam_hdr_destroy(header);
    sam_close(fp);

    return ClassifyReadType(platforms, readLengths);
}

bool ReadTypeUsesBaq(ReadType type) {
    return type != ReadType::Long;
}

int ResolveAdjustMQ(ReadType type, int requested) {
    if (requested != kAdjustMQByReadType) return requested;
    return type == ReadType::Long ? 0 : 40;
}
