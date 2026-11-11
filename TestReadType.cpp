/// Tests for read-type parsing, classification, detection and the pileup
/// defaults each read type selects. Detection tests write small SAM files to a
/// temporary directory, built inline so each test's input is visible here.

#include "ReadType.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

static int failures = 0;

#define EXPECT(cond)                                                                     \
    do {                                                                                 \
        if (!(cond)) {                                                                   \
            std::cerr << "  expectation failed: " #cond " (line " << __LINE__ << ")\n"; \
            ++failures;                                                                  \
        }                                                                                \
    } while (0)

static void Run(const char *name, void (*test)()) {
    const int before = failures;
    test();
    std::cout << (failures == before ? "PASS " : "FAIL ") << name << std::endl;
}

/// One SAM alignment line: `name` aligned at chr1:1 with an all-match CIGAR,
/// with its sequence written out or, if `storeSeq` is false, stored as '*'.
static std::string SamRecord(const std::string &name, int flag, int length, bool storeSeq = true) {
    return name + "\t" + std::to_string(flag) + "\tchr1\t1\t60\t" + std::to_string(length) + "M\t*\t0\t0\t" +
           (storeSeq ? std::string(length, 'A') : "*") + "\t*\n";
}

/// `count` records of `length` bp named with `prefix`.
static std::vector<std::string> SamRecords(const std::string &prefix, int count, int length) {
    std::vector<std::string> records;
    for (int i = 0; i < count; ++i) records.push_back(SamRecord(prefix + std::to_string(i), 0, length));
    return records;
}

/// Writes `header` and `records` to a new SAM file in a fresh temporary
/// directory and returns its path.
static std::string WriteSam(const std::string &header, const std::vector<std::string> &records) {
    const char *tmp = getenv("TMPDIR");
    std::string dirTemplate = std::string(tmp ? tmp : "/tmp") + "/vbid_readtype_XXXXXX";
    std::vector<char> dir(dirTemplate.begin(), dirTemplate.end());
    dir.push_back('\0');
    if (!mkdtemp(dir.data())) {
        std::cerr << "could not create a temporary directory\n";
        exit(EXIT_FAILURE);
    }
    const std::string path = std::string(dir.data()) + "/reads.sam";
    std::ofstream out(path);
    out << "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n" << header;
    for (const std::string &record : records) out << record;
    return path;
}

static void RemoveSam(const std::string &path) {
    remove(path.c_str());
    rmdir(path.substr(0, path.rfind('/')).c_str());
}

static ReadType Detect(const std::string &header, const std::vector<std::string> &records) {
    const std::string path = WriteSam(header, records);
    const ReadType type = DetectReadType(path.c_str(), nullptr).type;
    RemoveSam(path);
    return type;
}

static void parse_accepts_each_value_ignoring_case() {
    ReadType type = ReadType::Short;
    EXPECT(ParseReadType("auto", &type) && type == ReadType::Auto);
    EXPECT(ParseReadType("SHORT", &type) && type == ReadType::Short);
    EXPECT(ParseReadType("Long", &type) && type == ReadType::Long);
}

static void parse_rejects_unknown_value_and_leaves_output_unchanged() {
    ReadType type = ReadType::Long;
    EXPECT(!ParseReadType("ont", &type));
    EXPECT(!ParseReadType("", &type));
    EXPECT(type == ReadType::Long);
}

static void names_round_trip_through_parse() {
    for (ReadType type : {ReadType::Auto, ReadType::Short, ReadType::Long}) {
        ReadType parsed = ReadType::Auto;
        EXPECT(ParseReadType(ReadTypeName(type), &parsed) && parsed == type);
    }
}

static void short_reads_without_read_groups_are_short() {
    EXPECT(ClassifyReadType({}, std::vector<int32_t>(100, 150)).type == ReadType::Short);
}

static void median_at_threshold_is_long() {
    EXPECT(ClassifyReadType({}, {999, 1000, 1001}).type == ReadType::Long);
}

static void median_just_below_threshold_is_short() {
    EXPECT(ClassifyReadType({}, {998, 999, 50000}).type == ReadType::Short);
}

static void a_few_very_long_reads_do_not_make_short_reads_long() {
    std::vector<int32_t> lengths(90, 150);
    lengths.insert(lengths.end(), 10, 100000);
    EXPECT(ClassifyReadType({}, lengths).type == ReadType::Short);
}

static void ont_platform_is_long_even_for_short_reads() {
    EXPECT(ClassifyReadType({"ONT"}, std::vector<int32_t>(100, 150)).type == ReadType::Long);
}

static void platform_matching_ignores_case() {
    EXPECT(ClassifyReadType({"PacBio", "ont"}, {}).type == ReadType::Long);
}

static void illumina_platform_with_long_reads_follows_read_length() {
    EXPECT(ClassifyReadType({"ILLUMINA"}, std::vector<int32_t>(10, 20000)).type == ReadType::Long);
}

static void mixed_platforms_fall_back_to_read_length() {
    EXPECT(ClassifyReadType({"ONT", "ILLUMINA"}, std::vector<int32_t>(10, 150)).type == ReadType::Short);
}

static void a_read_group_without_platform_stops_platform_detection() {
    EXPECT(ClassifyReadType({"ONT", ""}, std::vector<int32_t>(10, 150)).type == ReadType::Short);
}

static void a_lone_read_group_without_platform_follows_read_length() {
    EXPECT(ClassifyReadType({""}, std::vector<int32_t>(10, 20000)).type == ReadType::Long);
}

static void no_platforms_and_no_reads_is_short() {
    EXPECT(ClassifyReadType({}, {}).type == ReadType::Short);
}

static void baq_is_used_only_for_short_reads() {
    EXPECT(ReadTypeUsesBaq(ReadType::Short));
    EXPECT(!ReadTypeUsesBaq(ReadType::Long));
}

static void adjust_mq_by_read_type_is_40_for_short_reads() {
    EXPECT(ResolveAdjustMQ(ReadType::Short, kAdjustMQByReadType) == 40);
}

static void adjust_mq_by_read_type_is_off_for_long_reads() {
    EXPECT(ResolveAdjustMQ(ReadType::Long, kAdjustMQByReadType) == 0);
}

static void an_explicit_adjust_mq_is_kept_for_long_reads() {
    EXPECT(ResolveAdjustMQ(ReadType::Long, 50) == 50);
}

static void an_explicit_adjust_mq_of_zero_is_kept_for_short_reads() {
    EXPECT(ResolveAdjustMQ(ReadType::Short, 0) == 0);
}

static void discard_near_indel_by_read_type_is_2_for_long_reads() {
    EXPECT(ResolveDiscardNearIndel(ReadType::Long, kDiscardNearIndelByReadType) == 2);
}

static void discard_near_indel_by_read_type_is_off_for_short_reads() {
    EXPECT(ResolveDiscardNearIndel(ReadType::Short, kDiscardNearIndelByReadType) == 0);
}

static void an_explicit_discard_near_indel_of_zero_is_kept_for_long_reads() {
    EXPECT(ResolveDiscardNearIndel(ReadType::Long, 0) == 0);
}

static void an_explicit_discard_near_indel_is_kept_for_short_reads() {
    EXPECT(ResolveDiscardNearIndel(ReadType::Short, 3) == 3);
}

static void each_read_counts_once_only_for_long_reads() {
    EXPECT(ReadTypeCountsEachReadOnce(ReadType::Long));
    EXPECT(!ReadTypeCountsEachReadOnce(ReadType::Short));
}

static void detects_long_reads_in_a_file_without_read_groups() {
    EXPECT(Detect("", {SamRecord("a", 0, 5000), SamRecord("b", 16, 8000), SamRecord("c", 0, 3000)}) ==
           ReadType::Long);
}

static void detects_short_reads_with_an_illumina_read_group() {
    EXPECT(Detect("@RG\tID:rg1\tSM:s1\tPL:ILLUMINA\n", SamRecords("r", 3, 150)) == ReadType::Short);
}

static void detects_long_reads_from_an_ont_read_group() {
    EXPECT(Detect("@RG\tID:rg1\tSM:s1\tPL:ONT\n", {SamRecord("a", 0, 200)}) == ReadType::Long);
}

static void detection_ignores_unmapped_secondary_and_supplementary_reads() {
    EXPECT(Detect("", {SamRecord("a", 4, 150), SamRecord("b", 256, 150), SamRecord("c", 2048, 150),
                       SamRecord("d", 0, 5000), SamRecord("e", 0, 6000)}) == ReadType::Long);
}

static void detection_ignores_duplicate_and_qc_failed_reads() {
    EXPECT(Detect("", {SamRecord("a", 1024, 150), SamRecord("b", 512, 150), SamRecord("c", 1024, 150),
                       SamRecord("d", 0, 5000)}) == ReadType::Long);
}

static void detection_measures_reads_stored_without_a_sequence() {
    EXPECT(Detect("", {SamRecord("a", 0, 5000, false), SamRecord("b", 0, 6000, false)}) == ReadType::Long);
}

static void detection_samples_only_the_first_reads() {
    std::vector<std::string> records = SamRecords("short", kReadTypeSampleSize, 150);
    const std::vector<std::string> longRecords = SamRecords("long", kReadTypeSampleSize, 20000);
    records.insert(records.end(), longRecords.begin(), longRecords.end());
    EXPECT(Detect("", records) == ReadType::Short);
}

static void a_file_with_only_a_header_is_short() {
    EXPECT(Detect("", {}) == ReadType::Short);
}

static void detection_fails_for_a_missing_file() {
    bool failed = false;
    try {
        DetectReadType("/nonexistent/reads.bam", nullptr);
    } catch (...) {
        failed = true;
    }
    EXPECT(failed);
}

int main() {
    Run("parse_accepts_each_value_ignoring_case", parse_accepts_each_value_ignoring_case);
    Run("parse_rejects_unknown_value_and_leaves_output_unchanged",
        parse_rejects_unknown_value_and_leaves_output_unchanged);
    Run("names_round_trip_through_parse", names_round_trip_through_parse);
    Run("short_reads_without_read_groups_are_short", short_reads_without_read_groups_are_short);
    Run("median_at_threshold_is_long", median_at_threshold_is_long);
    Run("median_just_below_threshold_is_short", median_just_below_threshold_is_short);
    Run("a_few_very_long_reads_do_not_make_short_reads_long", a_few_very_long_reads_do_not_make_short_reads_long);
    Run("ont_platform_is_long_even_for_short_reads", ont_platform_is_long_even_for_short_reads);
    Run("platform_matching_ignores_case", platform_matching_ignores_case);
    Run("illumina_platform_with_long_reads_follows_read_length", illumina_platform_with_long_reads_follows_read_length);
    Run("mixed_platforms_fall_back_to_read_length", mixed_platforms_fall_back_to_read_length);
    Run("a_read_group_without_platform_stops_platform_detection",
        a_read_group_without_platform_stops_platform_detection);
    Run("a_lone_read_group_without_platform_follows_read_length",
        a_lone_read_group_without_platform_follows_read_length);
    Run("no_platforms_and_no_reads_is_short", no_platforms_and_no_reads_is_short);
    Run("baq_is_used_only_for_short_reads", baq_is_used_only_for_short_reads);
    Run("adjust_mq_by_read_type_is_40_for_short_reads", adjust_mq_by_read_type_is_40_for_short_reads);
    Run("adjust_mq_by_read_type_is_off_for_long_reads", adjust_mq_by_read_type_is_off_for_long_reads);
    Run("an_explicit_adjust_mq_is_kept_for_long_reads", an_explicit_adjust_mq_is_kept_for_long_reads);
    Run("an_explicit_adjust_mq_of_zero_is_kept_for_short_reads",
        an_explicit_adjust_mq_of_zero_is_kept_for_short_reads);
    Run("discard_near_indel_by_read_type_is_2_for_long_reads", discard_near_indel_by_read_type_is_2_for_long_reads);
    Run("discard_near_indel_by_read_type_is_off_for_short_reads",
        discard_near_indel_by_read_type_is_off_for_short_reads);
    Run("an_explicit_discard_near_indel_of_zero_is_kept_for_long_reads",
        an_explicit_discard_near_indel_of_zero_is_kept_for_long_reads);
    Run("an_explicit_discard_near_indel_is_kept_for_short_reads",
        an_explicit_discard_near_indel_is_kept_for_short_reads);
    Run("each_read_counts_once_only_for_long_reads", each_read_counts_once_only_for_long_reads);
    Run("detects_long_reads_in_a_file_without_read_groups", detects_long_reads_in_a_file_without_read_groups);
    Run("detects_short_reads_with_an_illumina_read_group", detects_short_reads_with_an_illumina_read_group);
    Run("detects_long_reads_from_an_ont_read_group", detects_long_reads_from_an_ont_read_group);
    Run("detection_ignores_unmapped_secondary_and_supplementary_reads",
        detection_ignores_unmapped_secondary_and_supplementary_reads);
    Run("detection_ignores_duplicate_and_qc_failed_reads", detection_ignores_duplicate_and_qc_failed_reads);
    Run("detection_measures_reads_stored_without_a_sequence", detection_measures_reads_stored_without_a_sequence);
    Run("detection_samples_only_the_first_reads", detection_samples_only_the_first_reads);
    Run("a_file_with_only_a_header_is_short", a_file_with_only_a_header_is_short);
    Run("detection_fails_for_a_missing_file", detection_fails_for_a_missing_file);

    if (failures) {
        std::cerr << failures << " expectation(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "All read-type tests passed\n";
    return EXIT_SUCCESS;
}
