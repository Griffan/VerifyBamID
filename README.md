[![CI](https://github.com/Griffan/VerifyBamID/actions/workflows/ci.yml/badge.svg)](https://github.com/Griffan/VerifyBamID/actions/workflows/ci.yml)
[![GitHub Downloads](https://img.shields.io/github/downloads/Griffan/VerifyBamID/total.svg?style=flat)](https://github.com/Griffan/VerifyBamID/releases)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/verifybamid2/badges/downloads.svg)](https://anaconda.org/bioconda/verifybamid2)
[![DockerHub Pulls](https://img.shields.io/docker/pulls/griffan/verifybamid2)](https://hub.docker.com/r/griffan/verifybamid2)

# verifyBamID2

* Motivation: Detecting and estimating inter-sample DNA contamination became a crucial quality assessment step to ensure high quality sequence reads and reliable downstream analysis. Across many existing models, allele frequency usually is used to calculate prior genotype probability. Lack of accurate assignment of allele frequency could result in underestimation of contamination level. Hence we propose this ancestry-agnostic DNA contamination estimation method.

* Results: We applied our method to 1000 Genomes datasets by simulating contamination levels from 1% to 20% and comparing the contamination estimates obtain from different methods. When using pooled allele frequencies, as opposed to population-specific allele frequencies, we observed that the contamination levels are underestimated by 20%, 40%, 51%, and 73% for CEU, YRI, FIN, and CHS populations, respectively. Using our new method, the underestimation bias was reduced to 2-5%.

* Input Files: Aligned NGS sequence files(BAM or CRAM); Marker related files(SVD result on genotype matrix, provided in resorce directory)

* Please cite: 
Fan Zhang, Matthew Flickinger, Sarah A. Gagliano Taliun, InPSYght Psychiatric Genetics Consortium, Gonçalo R. Abecasis, Laura J. Scott, Steven A. McCaroll, Carlos N. Pato, Michael Boehnke, and Hyun Min Kang. 2020. “Ancestry-agnostic estimation of DNA sample contamination from sequence reads.” Genome Research. https://doi.org/10.1101/gr.246934.118

## Installation

The official project name is **verifyBamID2**. The command-line executable is named `VerifyBamID`.

### From source code
By default the build downloads [htslib](https://github.com/samtools/htslib) 1.25 and [libdeflate](https://github.com/ebiggers/libdeflate) 1.26 and links them into `VerifyBamID` statically, so it needs network access. htslib is always built with support for remote files (http, https, S3 and Google Cloud Storage) and for CRAM files compressed with lzma or bzip2. Its own dependencies are linked dynamically, so their development files must be installed, along with CMake 3.10 or later and a compiler that can build htslib 1.25 (clang, or GCC 8 or later; GCC 9 or later on 64-bit Arm Linux):

- Debian or Ubuntu: `apt-get install cmake g++ make libcurl4-openssl-dev libssl-dev zlib1g-dev libbz2-dev liblzma-dev`
- Fedora, RHEL or Amazon Linux: `dnf install cmake gcc gcc-c++ make libcurl-devel openssl-devel zlib-devel bzip2-devel xz-devel`
- macOS: the Xcode command-line tools provide everything except CMake, which Homebrew provides: `brew install cmake`

```
mkdir build
cd build
cmake ..
make
make test
```

The executable is written to `bin/VerifyBamID`. If the dependencies are installed somewhere CMake doesn't search, for example a conda environment, add `-DCMAKE_PREFIX_PATH=/that/prefix` to the `cmake` command.

#### Building without network access
To build the bundled htslib without network access, for example when packaging, download the [htslib 1.25](https://github.com/samtools/htslib/releases/download/1.25/htslib-1.25.tar.bz2) and [libdeflate 1.26](https://github.com/ebiggers/libdeflate/releases/download/v1.26/libdeflate-1.26.tar.gz) release tarballs beforehand and pass their absolute paths. Each may also be the tarball's extracted contents. A tarball must match the pinned checksum, and an extracted directory must hold the right version:

```
cmake -DHTSLIB_SOURCE=/path/to/htslib-1.25.tar.bz2 -DLIBDEFLATE_SOURCE=/path/to/libdeflate-1.26.tar.gz ..
```

#### Building against an installed htslib
To link an installed htslib 1.25 or later instead of building one, which then has to be present at run time too:

```
cmake -DUSE_SYSTEM_HTSLIB=ON ..
```

If htslib is not found, give its location, using `libhts.dylib` on macOS:

```
cmake -DUSE_SYSTEM_HTSLIB=ON -DHTS_INCLUDE_DIRS=/hts_absolute_path/include/ -DHTS_LIBRARIES=/hts_absolute_path/lib/libhts.so ..
```

Passing `HTS_INCLUDE_DIRS` or `HTS_LIBRARIES` the first time a build directory is configured turns on `USE_SYSTEM_HTSLIB` by itself.

### From bioconda
conda install -c bioconda verifybamid2 
  
## Guide for beginners

We have prepared example cmdline and necessary auxiliary resource files in the repository.
If you are unfamiliar with the arguments of verifyBamID2, you can verify the status of your BamFile by following examples listed below: 

### For BAM/CRAMs aligned to GRCh37 
(Note that GRCh37 assumes ***without chr prefix*** 1000 Genomes version - not UCSC - of human genome build 37, which is available at ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz ). Use the following commands for the BAM/CRAM files mapped to GRCh37.

```
$(VERIFY_BAM_ID_HOME)/bin/VerifyBamID \
  --UDPath $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b37.vcf.gz.dat.UD \
  --BedPath $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b37.vcf.gz.dat.bed \
  --MeanPath $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b37.vcf.gz.dat.mu \
  --Reference [/path/to/human_g1k_v37.fasta(.gz)] \
  --BamFile [/path/to/bam/or/cram/file]
```
or simply:
```
$(VERIFY_BAM_ID_HOME)/bin/VerifyBamID \
  --SVDPrefix $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b37.vcf.gz.dat \
  --Reference [/path/to/human_g1k_v37.fasta(.gz)] \
  --BamFile [/path/to/bam/or/cram/file]
```
### For BAM/CRAMs aligned to GRCh38
(Note that GRCh38 assumes ***with chr prefix*** 1000 Genomes version of human genome build 38, which is available at ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome ). Use the following commands for the BAM/CRAM files mapped to GRCh38.

```
$(VERIFY_BAM_ID_HOME)/bin/VerifyBamID \
  --UDPath $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b38.vcf.gz.dat.UD \
  --BedPath $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b38.vcf.gz.dat.bed \
  --MeanPath $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b38.vcf.gz.dat.mu \
  --Reference [/path/to//GRCh38_full_analysis_set_plus_decoy_hla.fa] \
  --BamFile [/path/to/bam/or/cram/file] 
```
or simply:
```
$(VERIFY_BAM_ID_HOME)/bin/VerifyBamID \
  --SVDPrefix $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b38.vcf.gz.dat \
  --Reference [/path/to//GRCh38_full_analysis_set_plus_decoy_hla.fa] \
  --BamFile [/path/to/bam/or/cram/file] 
```

### Resource files are ready
(for both 1000 Genome Project(1000g) dataset and Human Genome Diversity Project(hgdp) dataset are available)

You can directly use reference panel information by using our pre-calculated resource files in ``$(VERIFY_BAM_ID_HOME)/resource/`` directory.

## Usage
For regular estimation:
```
$(VERIFY_BAM_ID_HOME)/bin/VerifyBamID --BamFile [/path/to/bam/or/cram/file] --UDPath [/path/to/UD/file] --BedPath [/path/to/bed/file] --MeanPath [/path/to/mu/file] --Reference [/path/to/fasta/file]
```
or simply:
```
$(VERIFY_BAM_ID_HOME)/bin/VerifyBamID --BamFile [/path/to/bam/or/cram/file] --SVDPrefix [/path/to/SVDPrefix/file] --Reference [/path/to/fasta/file]
```

```
--SVDPrefix      [String] SVD related files prefix(normally shared by .UD, .mu and .bed files)[Required]
--BamFile        [String] Bam or Cram file for the sample[Required]
--Reference      [String] reference file[Required]
--Seed           [INT] Random number seed(default:12345)
--NumPC          [INT] Number of Principal Components used in estimation
--NumThread      [Int] Set number of threads in likelihood calculation[default:4]
--FixPC          [String] Specify known PC coordinates for the sample[format PC1:PC2:PC3...]
--FixAlpha       [Double] Specify known contamination level
--WithinAncestry [Bool] Enabling withinAncestry assume target sample and contamination source are from the same populations,[default:betweenAncestry] otherwise")
--KnownAF        [String] A Bed file that provide known allele frequency for each marker, similar behaviour with VerifyBamID 1.0
--Epsilon        [Double] Minimization procedure convergence threshold, usually a trade-off bettween accuracy and running time[default:1e-10]
--OutputPileup   [Bool] If output temp pileup file
--Verbose        [Bool] If print the progress of the method on the screen
/*To construct SVDPrefix auxillary files*/
--RefVCF         [String] Reference panel VCF with genotype information, for generation of .UD .mu .bed files[Optional]
/*Pileup information related*/
--ReadType       [String] auto, short or long. Selects pileup defaults for the reads; see "Long reads" below. [default:auto]
--Methylation    [Bool] Directional bisulfite (WGBS) / enzymatic methyl-seq (EM-seq) mode. Uses only observations that unambiguously identify an allele on each read's bisulfite-conversion strand, and disables MAPQ adjustment (--adjust-MQ) and BAQ (both penalize the expected C-to-T/G-to-A conversions). See the "Methylation sequencing" section below. [Optional]
--no-orphans     [Bool] Skip anomalous read pairs in variant calling. Anomolous read pairs are those marked in the FLAG field as paired in sequencing but without the properly-paired flag set. [Optional]
--ignoreRG       [Bool] Ignore read-group tags and treat all reads as one sample group regardless of BAM read-group headers. Reads are grouped by the input file name. [Optional]
--adjust-MQ      [Int] Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. A zero value disables this functionality; if enabled, the recommended value for BWA is 50. -1 means 40 for short reads and 0 for long reads. [default:-1]
--max-depth      [Int] Setting this limit reduces the amount of memory and time needed to process regions with very high coverage. Passing zero for this option sets it to the highest possible value, effectively removing the depth limit. [8000]
--incl-flags     [Int] Required flags: skip reads with mask bits unset [null]
--excl-flags     [Int] Filter flags: skip reads with mask bits set [UNMAP,SECONDARY,QCFAIL,DUP]
--DiscardNearIndel [Int] Discard bases this many or fewer aligned bases from an insertion or deletion in their read. 0 disables; -1 means 2 for long reads and 0 for short reads. [default:-1]
/*Below are deprecated but still available*/
--UDPath         [String] .UD matrix file from SVD result of genotype matrix[Required]
--MeanPath       [String] .mu matrix file of genotype matrix[Required]
--BedPath        [String] .Bed file for markers used in this analysis,format(chr\tpos-1\tpos\trefAllele\taltAllele)[Required]
```

## Output Format

```
Estimation from OptimizeHeter:
Contaminating Sample PC1:-0.623602      PC2:0.57292
Intended Sample  PC1:-0.036304  PC2:0.0200112
Alpha:0.0013662
```
First line: Which model used

Second line: PC coordinates of Contaminating Sample

Third line: PC coordinates of Intended Sample(the sample in which you are interested)

Fourth line: Estimated contamination level 


Also, upon the completion of each run, you may find two files ending with suffix:
* “.selfSM” which shares the same format as VB1(https://genome.sph.umich.edu/wiki/VerifyBamID), and the key information FREEMIX indicates the estimated contamination level.
* “.Ancestry” which contains the PC coordinates for both intended sample and contaminating sample, with each row being one PC.

## Long reads

VerifyBamID works with long reads (Oxford Nanopore, PacBio) as well as short reads. `--ReadType` selects suitable pileup settings, and by default (`auto`) it detects the read type from the BAM/CRAM:

* The reads are **long** if every `@RG` header line has `PL:ONT` or `PL:PACBIO`, or if the median length of the first 500 primary mapped reads (excluding duplicates and QC failures) is at least 1000 bp. Otherwise they are **short**. The decision and its reason are logged. `--ReadType` and `--DiscardNearIndel` have no effect with `--PileupFile`.
* For **long** reads, BAQ is disabled and `--adjust-MQ` defaults to 0. BAQ's cost grows with read length times the read's net insertion/deletion imbalance, which for long reads can mean hours of computation and tens of gigabytes of memory per read. The `--adjust-MQ` mismatch heuristic was designed for short reads.
* Also for **long** reads, bases within 2 aligned bases of an insertion or deletion in their read are discarded (`--DiscardNearIndel` defaults to 2). When an indel error sits next to a SNP, the aligner often places it on the far side and reads the SNP as the reference base, which looks like contamination. On clean public ONT and PacBio samples this cuts the reported contamination about fivefold, at a cost of 0.5–1.3% of the bases.
* Also for **long** reads, each read counts at most once at a marker. An ONT fold-back read can align twice over the same position, as a primary alignment and a supplementary alignment on the other strand, and both copies observe the same molecule. If the copies show different alleles (a deletion counts as one), all of them are discarded, including copies that the base-quality or near-indel filters would have discarded anyway. Otherwise the copy with the highest base quality that passes those filters is kept.
* For **short** reads, BAQ is on, `--adjust-MQ` defaults to 40, no bases are discarded near indels (BAQ already lowers their qualities), and overlapping mates are handled as in `samtools mpileup`.

An explicit `--adjust-MQ` or `--DiscardNearIndel` value is always used as given. `--Methylation` still disables both BAQ and MAPQ adjustment; with short reads that leaves bases near indels unprotected, so consider `--DiscardNearIndel 2` there. Counting each read once has no separate switch: `--ReadType short` turns it off along with the other long-read settings. Use `--ReadType short` or `--ReadType long` to override detection.

**Choose a panel whose markers are further apart than your reads are long.** Each read samples a single molecule from either the intended sample or the contaminant, so it is one piece of evidence. A read that covers several markers is counted at each of them, so the likelihood overstates the evidence, and where markers cluster, reads over the clusters count for more than others. Aim for each read to cover, on average, no more than one marker. The bundled 10k panels (for example `resource/1000g.phase3.10k.b38.vcf.gz.dat` or `resource/hgdp.10k.b38.vcf.gz.dat`) have markers about 270 kb apart on average (median about 150 kb) and suit long reads. The 100k panels have markers about 27 kb apart (median about 14 kb), so long reads often cover several markers. When building a panel with `--RefVCF`, thin the sites to a similar spacing.

## Methylation sequencing (bisulfite / EM-seq)

In whole-genome bisulfite sequencing (WGBS) and enzymatic methyl-seq (EM-seq), unmethylated cytosines are read as `T` (and, on reads from the opposite strand, unmethylated `G` positions read as `A`). A regular pileup therefore disagrees with the reference at a large fraction of positions for reasons unrelated to genotype, which would inflate the contamination estimate. The `--Methylation` flag enables a mode that handles this correctly, reusing the same SVD reference panels.

```
$(VERIFY_BAM_ID_HOME)/bin/VerifyBamID \
  --Methylation \
  --SVDPrefix $(VERIFY_BAM_ID_HOME)/resource/1000g.phase3.100k.b38.vcf.gz.dat \
  --Reference [/path/to/GRCh38.fa] \
  --BamFile [/path/to/methyl-seq.bam]
```

### What it does

At each marker, `--Methylation` uses only the observations that still identify an allele unambiguously on the read's conversion strand, and discards the rest:

* **A/T** markers are used from both strands (conversion never touches A or T).
* **A/C and G/T** markers are used from the one strand on which conversion cannot mimic the other allele (half depth).
* **A/G and C/T** (transitions) are used only from the single strand that stays informative — A/G from C-to-T-converted reads, C/T from G-to-A-converted reads (half depth).
* **C/G** markers are **not used**: both strands are ambiguous (C is affected on the C-to-T strand, G on the G-to-A strand), so no observations survive. This affects 2% of the 1000 Genomes panel and 0% of HGDP (which excludes palindromic SNPs). When building a custom marker panel with `--RefVCF`, C/G markers can be excluded up front.
* Observations that a conversion could make ambiguous are dropped rather than guessed at.

Because ambiguous observations are dropped rather than disambiguated, **the estimate does not depend on the methylation level or on the bisulfite/enzymatic conversion efficiency**, and **CpG context is irrelevant** — no methylation-rate or conversion-rate parameter enters the model. The cost is reduced effective depth (roughly half the observations at transition and most transversion markers), so somewhat higher coverage is recommended than for a standard run.

`--Methylation` also changes read handling to suit converted reads:

* **MAPQ adjustment (`--adjust-MQ`) and BAQ are disabled.** Both score reads against the *unconverted* reference, so an EM-seq/WGBS read's conversions look like dense mismatches; left on, MAPQ capping rejects the great majority of reads outright.
* **Supplementary alignments are excluded**, and a read is used only when its conversion strand can be determined reliably — from an aligner methylation tag (`YD`, `ZS`, or `XG`) when present, otherwise from a proper-pair FLAG (paired reads must be properly paired; single-end reads use their own strand). A read whose strand cannot be determined is skipped, because a misassigned strand would read converted bases as alternate alleles and inflate the estimate.

### Scope and limitations

* **Directional libraries only.** This covers EM-seq and standard directional (Lister-style) WGBS. Non-directional and PBAT libraries are supported only where the aligner recorded a strand tag; with FLAG-only inference their strand assignment is wrong, so the estimate would be unreliable. Strand is taken from an aligner methylation tag (bismark `XG`, bwa-meth/BISCUIT `YD`, BSMAP `ZS`) when present, otherwise inferred from read-pair orientation. The mode has been validated end-to-end on directional EM-seq (FLAG-inferred strand); tag-based strand assignment for the other aligners is covered by unit tests but not yet exercised on their real output.
* **`--PileupFile` input.** A pileup text file does not carry read-pair or tag information, so the conversion strand cannot be recovered from it. `--Methylation --PileupFile` is allowed but only for input that has *already* been filtered to methylation-safe observations (for example, a pileup produced by an earlier `--Methylation --BamFile ... --OutputPileup` run); a warning is printed. Feeding an unfiltered pileup will inflate the estimate.
* **TAPS+ and Illumina 5-base chemistries.** Methods such as TAPS+ and Illumina's 5-base prep convert **methylated** C to T rather than unmethylated C. Because the substitution direction is the same (C→T on the top strand, G→A on the bottom), `--Methylation` mode is applicable and will produce a correct estimate. However, these chemistries convert far fewer positions per read (predominantly CpG-context cytosines) than bisulfite/EM-seq, so `--Methylation` mode is likely more conservative than necessary — it drops observations at all C- or G-containing markers on the affected strand, when in practice only CpG-context markers are at significant risk. For datasets with sufficient coverage `--Methylation` will likely give the cleanest estimate, but the mode has not been broadly tested on TAPS+ or 5-base data. A future version may add a dedicated mode for methylated-C→T chemistries with CpG-context-aware filtering and tuned MAPQ/BAQ thresholds.

## Running from docker

To run VB2 from docker hub, you can simply:
```
docker pull griffan/verifybamid2

```
and run the image by mounting local directory:
```
docker run -v  /localPath/localDirectory:/VerifyBamID/localDirectory  griffan/verifybamid2:v1.0.6 VerifyBamID --SVDPrefix /VerifyBamID/resource/1000g.phase3.10k.b37.vcf.gz.dat --BamFile  /VerifyBamID/localDirectory/targetSample.bam --Reference  /VerifyBamID/localDirectory/hs37d5.fa 
```
Make sure the "targetSample.bam" and "hs37d5.fa" are located inside "/localPath/localDirectory"

## Generating your own resource files.

For producing customized resource files to be used as the input argument of verifyBamID2 generation, you need to start with a VCF file and FASTA formatted reference files. Please refer to the example below.
```
VerifyBamID --RefVCF ReferencePanel.vcf.gz --Reference ./resource/test/chr20.fa.gz

```
Note that for the ReferencePanel.vcf.gz, VerifyBamID expects each marker:
* bi-allelic 
* SNP only 
* with "PASS" filter label 
* with valid GT,GL or PL field
* genotype missing rate less than 20%
* at least 5000 markers

In the example above, the expected output file names will be ``ReferencePanel.vcf.gz.UD, ReferencePanel.vcf.gz.mu, ReferencePanel.vcf.gz.bed``

## Need more markers

If you want to use more than 1M markers, the memory consumption of verifyBamID2 could be large. One way to circumvent this situation is to use 'a truncated SVD' technique, which is described in this link:
http://bwlewis.github.io/1000_genomes_examples/PCA_whole_genome.html

## Generating PC plot

After each run, you will get the contamination Alpha estimation, as well as ancestry PC coordinates for both intended sample and contaminating sample.

You may want to visualize these information, in that case, the PC coordinates files(ending with .V) in ``$(VERIFY_BAM_ID_HOME)/resource/`` might help you by
providing background PC points of 1000 Genomes Project samples(e.g. 1000g.100k.b38.vcf.gz.dat.V) or of Human Genome Diversity Project samples(e.g. hgdp.100k.b38.vcf.gz.dat.V)

We also provide script to generate PC plot with customized dataset as background points, for example:
```
sh $(VERIFY_BAM_ID_HOME)/bin/run.plot.sh -i ./resource/test/hapmap_3.3.b37.dat.V -o ./resource/test/hapmap -r 1000g -g grey
```
You may run ``sh $(VERIFY_BAM_ID_HOME)/bin/run.plot.sh -h`` for further help.

## Debugging Mode

If you encounter abnormal errors, e.g. "Segmentation fault", you can try to recompile the build under debugging mode:
```
cmake .. -DCMAKE_BUILD_TYPE=Debug
```
and then rerun your command line to locate the backtrace message, and then post it to issues page, for example:
```
Stack trace (most recent call last):
#6    Object "VerifyBamID", at 0x10ad9d822, in main + 466
#5    Object "VerifyBamID", at 0x10ad9b53d, in execute(int, char**) + 8509
#4    Object "VerifyBamID", at 0x10adc70c7, in ContaminationEstimator::OptimizeLLK(std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> > const&) + 903
#3    Object "VerifyBamID", at 0x10adc9408, in ContaminationEstimator::OptimizeHeter(AmoebaMinimizer&) + 1560
#2    Object "libsystem_platform.dylib", at 0x7ff8182f0dfc, in _sigtramp + 28
#1    Object "VerifyBamID", at 0x10ad9e3ed, in backward::SignalHandling::sig_handler(int, __siginfo*, void*) + 13
#0    Object "VerifyBamID", at 0x10ad9e456, in backward::SignalHandling::handleSignal(int, __siginfo*, void*) + 70
```
## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## Bug Report

fanzhang@umich.edu

## License

This project is licensed under the terms of the MIT license.
