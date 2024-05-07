[![Build Status](https://travis-ci.org/Griffan/VerifyBamID.png?branch=master)](https://travis-ci.org/Griffan/VerifyBamID)
[![GitHub Downloads](https://img.shields.io/github/downloads/Griffan/VerifyBamID/total.svg?style=flat)](https://github.com/Griffan/VerifyBamID/releases)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/verifybamid2/badges/downloads.svg)](https://anaconda.org/bioconda/verifybamid2)
[![DockerHub Pulls](https://img.shields.io/docker/pulls/griffan/verifybamid2)](https://hub.docker.com/r/griffan/verifybamid2)

# VerifyBamID2

* Motivation: Detecting and estimating inter-sample DNA contamination became a crucial quality assessment step to ensure high quality sequence reads and reliable downstream analysis. Across many existing models, allele frequency usually is used to calculate prior genotype probability. Lack of accurate assignment of allele frequency could result in underestimation of contamination level. Hence we propose this ancestry-agnostic DNA contamination estimation method.

* Results: We applied our method to 1000 Genomes datasets by simulating contamination levels from 1% to 20% and comparing the contamination estimates obtain from different methods. When using pooled allele frequencies, as opposed to population-specific allele frequencies, we observed that the contamination levels are underestimated by 20%, 40%, 51%, and 73% for CEU, YRI, FIN, and CHS populations, respectively. Using our new method, the underestimation bias was reduced to 2-5%.

* Input Files: Aligned NGS sequence files(BAM or CRAM); Marker related files(SVD result on genotype matrix, provided in resorce directory)

* Please cite: 
Fan Zhang, Matthew Flickinger, Sarah A. Gagliano Taliun, InPSYght Psychiatric Genetics Consortium, Gonçalo R. Abecasis, Laura J. Scott, Steven A. McCaroll, Carlos N. Pato, Michael Boehnke, and Hyun Min Kang. 2020. “Ancestry-agnostic estimation of DNA sample contamination from sequence reads.” Genome Research. https://doi.org/10.1101/gr.246934.118

## Installation

### From source code
  - mkdir build
  - cd build
  - cmake ..

```
In case any required libraries is missing, you may specify customized installing path by replacing "cmake .." with:

For libhts:
  - cmake -DHTS_INCLUDE_DIRS=/hts_absolute_path/include/  -DHTS_LIBRARIES=/hts_absolute_path/lib/libhts.a ..

For bzip2:
  - cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

For lzma:
  - cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..
```

  - make
  - make test
  
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
--no-orphans     [Bool] Skip anomalous read pairs in variant calling. Anomolous read pairs are those marked in the FLAG field as paired in sequencing but without the properly-paired flag set. [Optional]
--adjust-MQ      [Int] Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. A zero value disables this functionality; if enabled, the recommended value for BWA is 50. [default:40]
--max-depth      [Int] Setting this limit reduces the amount of memory and time needed to process regions with very high coverage. Passing zero for this option sets it to the highest possible value, effectively removing the depth limit. [8000]
--incl-flags     [Int] Required flags: skip reads with mask bits unset [null]
--excl-flags     [Int] Filter flags: skip reads with mask bits set [UNMAP,SECONDARY,QCFAIL,DUP]
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
