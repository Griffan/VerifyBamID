# Changelog

All notable changes to verifyBamID2 are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Dates are those of the GitHub releases.

## [2.1.0] (pending)

### Added

- `--ReadType auto|short|long` chooses pileup settings suited to the reads ([#99]). `auto`, the default, chooses long when every `@RG` has `PL:ONT` or `PL:PACBIO`, or when the median length of the first 500 primary mapped reads is at least 1000 bp. Long reads get BAQ off, `--adjust-MQ 0` and `--DiscardNearIndel 2`, and each read counts once per marker. Short reads keep the previous settings.
- `--DiscardNearIndel N` discards bases within N aligned bases of an insertion or deletion in their read ([#99]). It defaults to 2 for long reads and 0 (off) for short reads.
- `--Methylation` estimates contamination in bisulfite and EM-seq libraries, using only observations that C→T and G→A conversion cannot mimic ([#90]). It supports directional libraries.
- `--ignoreRG` treats all reads as one sample, whatever their read groups ([#92]).

### Changed

- BAM and CRAM input is read in one pass over all markers instead of one index query per marker, and CRAM decodes only the fields the pileup uses ([#99]). Long-read genomes that never finished now take minutes: a 27× ONT genome on the 100k panel takes about 7 minutes and 1.4 GB. A 30× Illumina CRAM went from 8.2 to 4.6 minutes, with identical output.
- htslib 1.25 or later is required ([#99]). The default build now downloads htslib 1.25 and libdeflate 1.26, checks their SHA-256, and links them statically ([#100]). `-DUSE_SYSTEM_HTSLIB=ON` links an installed htslib instead, and `-DHTSLIB_SOURCE` and `-DLIBDEFLATE_SOURCE` build from local copies without network access.
- CMake 3.10 or later is required ([#100]).
- The Docker image is built on Ubuntu 26.04 with the bundled htslib, and a `REPO` build argument selects the repository to build from ([#87], [#100]).

### Fixed

- `--max-depth` applied only at the first marker; later markers used htslib's default of 8000 ([#97]). `--max-depth 0` now means no limit, as documented, rather than about one read per position, and negative values are rejected.
- Markers deep enough for their likelihood to underflow were silently left out of the estimate ([#94], [#98]).
- At a marker where a read had a deletion, the bases of later reads were paired with the wrong base qualities ([#90]).
- `--incl-flags` overwrote pileup options instead of setting the required-flag mask, so `--incl-flags 0` turned off overlapping-mate correction and `--incl-flags 2` did not require proper pairs ([#95]).
- `--RefVCF` no longer aborts on VCF lines whose sample columns omit trailing FORMAT fields, as GATK CombineGVCFs writes them; the omitted fields are read as missing ([#93]).

## [2.0.3] - 2026-06-29

### Added

- `--GramSVD` computes the reference-panel SVD from the Gram matrix, using far less memory on panels with many variants ([#83]). On a panel of about 950 samples and 2.5 million variants, it cut the SVD from hours to about 20 minutes on one CPU.
- The log reports how long each phase of estimation takes ([#82]).

### Changed

- Likelihood optimisation is about 20 times faster, with bit-identical results ([#82]).

### Fixed

- `--RefVCF` read every sample's genotype from the same FORMAT field, so a sample whose PL was `.` counted as missing even when it had a GT. Each sample now falls back from PL to GL to GT ([#76]).
- `--RefVCF` could not read genotypes given only as GT, such as `0/1` ([#76]).
- `--RefVCF` gave a sample with no readable genotype the previous sample's genotype ([#76]).

## [2.0.2] - 2026-04-11

### Added

- `--IncludeChr` sets the chromosomes `--RefVCF` uses, for non-human panels ([#75]). The default is still the human autosomes.
- `--SkipMinSampleCountCheck` lets `--RefVCF` use a panel with fewer than 1000 samples, with a warning ([#77]).
- `--NumSVDPCs` sets how many principal components `--RefVCF` writes (default 10; 0 writes all), and the log reports the variance each component explains ([#78]).
- `--RefVCF` logs its progress and the number of markers on each chromosome ([#79]).
- Debug builds (`-DCMAKE_BUILD_TYPE=Debug`) print a stack trace on a crash ([#65]).

### Changed

- `--RefVCF` skips markers that are not PASS, not biallelic or not SNPs, with a warning for each. A duplicated marker is an error.
- `--NumPC` is checked against the number of PCs in the `--SVDPrefix` files, instead of a fixed limit of 4.
- `--max-depth` values below 8000 are no longer raised to 8000 divided by the number of samples.
- The b38 exome resource files name chromosomes with a `chr` prefix, matching GRCh38.
- CMake 3.5 or later is required.

### Fixed

- `--PileupFile` input containing insertions or deletions paired bases with the wrong qualities, giving wrong estimates or a crash ([#59]).
- `--RefVCF` misread samples whose PL or GL field was `.`.
- The bundled Eigen 3.3.0 failed to compile with recent Clang ([#74]).
- The Docker image failed to build htslib ([#31]).

## [2.0.1] - 2021-01-13

### Added

- `--min-BQ`, `--min-MQ`, `--adjust-MQ`, `--max-depth`, `--incl-flags`, `--excl-flags` and `--no-orphans` set the pileup options.
- Exome marker resources for 1000 Genomes phase 3, on b37 and b38.

### Changed

- The version number jumps to 2.0.1 to distinguish verifyBamID2 from VerifyBamID 1.
- `--RefVCF` writes only the first 10 dimensions of the `.UD` and `.V` matrices.
- The format of the PCs in `.Ancestry` files and of the verbose output.

### Fixed

- Compilation on macOS Catalina.
- `hgdp.10k.b38.vcf.gz.dat.bed` had a redundant chromosome prefix.

## [1.0.6] - 2019-07-30

### Added

- A Dockerfile ([#4]).
- The version is printed at start-up.

### Changed

- htslib must be installed before building; the build no longer downloads it. `HTS_INCLUDE_DIRS` and `HTS_LIBRARIES` point CMake at it.

### Removed

- The 1000 Genomes phase 1 resource files.

### Fixed

- A bug when the sanity check was disabled with `--DisableSanityCheck`.
- The `hgdp.100k.b37` resource files, and the dimensions of the HGDP `.V` files.

## [1.0.5] - 2018-06-13

### Added

- The log reports the number of markers in the reference matrix.
- `--RefVCF` warns about markers outside the autosomes.

### Fixed

- Parsing of the sample name.
- The reported total number of input markers.

## [1.0.4] - 2018-04-22

### Fixed

- The sanity check failed on `--PileupFile` input.

## [1.0.3] - 2018-04-20

### Added

- `--SVDPrefix` names the `.UD`, `.mu` and `.bed` files by their shared prefix.
- 1000 Genomes phase 3 SVD resource files.

### Changed

- `--NumPC` above 4 is an error, since the bundled resource files have 4 PCs.
- Marker filtering conditions were adjusted.

### Fixed

- Bases were translated incorrectly from `--PileupFile` input.
- Failures writing output files are detected.

## [1.0.2] - 2018-01-14

### Added

- `--PileupFile` reads a pre-computed pileup instead of a BAM or CRAM.

### Changed

- Orphan reads are used by default.
- The result format is backward compatible with earlier releases.

### Fixed

- The format of `--OutputPileup` files.

## [1.0.1] - 2017-06-15

### Added

- Output files compatible with VerifyBamID 1, including `.selfSM`.
- `--Epsilon` sets the convergence threshold, trading speed against accuracy.

### Changed

- The output format.

## [1.0.0] - 2017-04-01

The first stable release.

[2.1.0]: https://github.com/Griffan/VerifyBamID/compare/v2.0.3...HEAD
[2.0.3]: https://github.com/Griffan/VerifyBamID/compare/V2.0.2...v2.0.3
[2.0.2]: https://github.com/Griffan/VerifyBamID/compare/2.0.1...V2.0.2
[2.0.1]: https://github.com/Griffan/VerifyBamID/compare/1.0.6...2.0.1
[1.0.6]: https://github.com/Griffan/VerifyBamID/compare/1.0.5...1.0.6
[1.0.5]: https://github.com/Griffan/VerifyBamID/compare/1.0.4...1.0.5
[1.0.4]: https://github.com/Griffan/VerifyBamID/compare/1.0.3...1.0.4
[1.0.3]: https://github.com/Griffan/VerifyBamID/compare/1.0.2...1.0.3
[1.0.2]: https://github.com/Griffan/VerifyBamID/compare/1.0.1...1.0.2
[1.0.1]: https://github.com/Griffan/VerifyBamID/compare/v1.0.0...1.0.1
[1.0.0]: https://github.com/Griffan/VerifyBamID/releases/tag/v1.0.0
[#4]: https://github.com/Griffan/VerifyBamID/pull/4
[#31]: https://github.com/Griffan/VerifyBamID/pull/31
[#59]: https://github.com/Griffan/VerifyBamID/pull/59
[#65]: https://github.com/Griffan/VerifyBamID/pull/65
[#74]: https://github.com/Griffan/VerifyBamID/pull/74
[#75]: https://github.com/Griffan/VerifyBamID/pull/75
[#76]: https://github.com/Griffan/VerifyBamID/pull/76
[#77]: https://github.com/Griffan/VerifyBamID/pull/77
[#78]: https://github.com/Griffan/VerifyBamID/pull/78
[#79]: https://github.com/Griffan/VerifyBamID/pull/79
[#82]: https://github.com/Griffan/VerifyBamID/pull/82
[#83]: https://github.com/Griffan/VerifyBamID/pull/83
[#87]: https://github.com/Griffan/VerifyBamID/pull/87
[#90]: https://github.com/Griffan/VerifyBamID/pull/90
[#92]: https://github.com/Griffan/VerifyBamID/pull/92
[#93]: https://github.com/Griffan/VerifyBamID/pull/93
[#94]: https://github.com/Griffan/VerifyBamID/pull/94
[#95]: https://github.com/Griffan/VerifyBamID/pull/95
[#97]: https://github.com/Griffan/VerifyBamID/pull/97
[#98]: https://github.com/Griffan/VerifyBamID/pull/98
[#99]: https://github.com/Griffan/VerifyBamID/pull/99
[#100]: https://github.com/Griffan/VerifyBamID/pull/100
