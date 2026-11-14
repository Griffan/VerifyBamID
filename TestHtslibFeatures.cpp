/// Checks that the htslib VerifyBamID is linked against was built with each
/// feature named on the command line, as reported by hts_features(). Exits 0
/// if all are present, 1 if any is missing, and 2 for an unknown name.

#include <htslib/hts.h>

#include <cstring>
#include <iostream>

struct Feature {
    const char *name;
    unsigned int bit;
};

static const Feature kFeatures[] = {
    {"lzma", HTS_FEATURE_LZMA},       {"bzip2", HTS_FEATURE_BZIP2}, {"libcurl", HTS_FEATURE_LIBCURL},
    {"s3", HTS_FEATURE_S3},           {"gcs", HTS_FEATURE_GCS},     {"libdeflate", HTS_FEATURE_LIBDEFLATE},
};

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " feature..." << std::endl;
        return 2;
    }
    const unsigned int available = hts_features();
    int status = 0;
    for (int i = 1; i < argc; ++i) {
        const Feature *feature = nullptr;
        for (const Feature &candidate : kFeatures) {
            if (std::strcmp(argv[i], candidate.name) == 0) feature = &candidate;
        }
        if (feature == nullptr) {
            std::cerr << "unknown htslib feature: " << argv[i] << std::endl;
            return 2;
        }
        const bool present = (available & feature->bit) != 0;
        std::cout << (present ? "present " : "MISSING ") << feature->name << std::endl;
        if (!present) status = 1;
    }
    std::cout << "htslib " << hts_version() << ": " << hts_feature_string() << std::endl;
    return status;
}
