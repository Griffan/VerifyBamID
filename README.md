# VerifyBamID

* Motivation: Detecting and estimating inter-sample DNA contamination became a crucial quality assessment step to ensure high quality sequence reads and reliable downstream analysis. Across many existing models, allele frequency usually is used to calculate prior genotype probability. Lack of accurate assignment of allele frequency could result in underestimation of contamination level.

* Results: We applied our method to 1000 Genomes datasets by simulating contamination levels from 1% to 20% and comparing the contamination estimates obtain from different methods. When using pooled allele frequencies, as opposed to population-specific allele frequencies, we observed that the contamination levels are underestimated by 20%, 40%, 51%, and 73% for CEU, YRI, FIN, and CHS populations, respectively. Using our new method, the underestimation bias was reduced to 2-5%.

* Input Files: Aligned NGS sequence files(Bam or Cram); Marker related files(SVD result on genotype matrix, provided in resorce directory)


## Installation

  - mkdir build
  - cd build
  - cmake ..
  - make
  - make test
  
## Guide for beginners

If you are unfamiliar with detailed arguments of verifyBamID2, use the resources provided in the repository as input files. More specifically try to run the following commands

### For BAM/CRAMs aligned to GRCh37 
(Note that GRCh37 assumes ***without chr prefix*** 1000 Genomes version - not UCSC - of human genome build 37, which is available at ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz ). Use the following commands for the BAM/CRAM files mapped to GRCh37.

```
$(VERIFY_BAM_ID_HOME)/bin/VerifyBamID --UDPath $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b37.vcf.gz.dat.UD \
  --BamFile [/path/to/bam/or/cram/file] --BedPath $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b37.vcf.gz.dat.bed \
  --MeanPath $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b37.vcf.gz.dat.mu \
  --Reference [/path/to/human_g1k_v37.fasta(.gz)]
```

### For BAM/CRAMs aligned to GRCh38
(Note that GRCh38 assumes ***with chr prefix*** 1000 Genomes version of human genome build 38, which is available at ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome ). Use the following commands for the BAM/CRAM files mapped to GRCh37.

```
$(VERIFY_BAM_ID_HOME)/bin/VerifyBamID --UDPath $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b38.vcf.gz.dat.UD \
   --BamFile [/path/to/bam/or/cram/file] --BedPath $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b38.vcf.gz.dat.bed \
   --MeanPath $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b38.vcf.gz.dat.mu \
   --Reference [/path/to//GRCh38_full_analysis_set_plus_decoy_hla.fa]
```

## Usage
For regular estimation:
```
$(VERIFY_BAM_ID_HOME)/binVerifyBamID --UDPath [/path/to/UD/file] --BamFile [/path/to/bam/or/cram/file] --BedPath [/path/to/bed/file] --MeanPath [/path/to/mu/file] --Reference [/path/to/fasta/file]
```
```
--UDPath    [String] .UD matrix file from SVD result of genotype matrix[Required]
--MeanPath  [String] .Mean matrix file of genotype matrix[Required]
--BedPath   [String] .Bed file for markers used in this analysis,format(chr\tpos-1\tpos\trefAllele\taltAllele)[Required]
--RefVCF    [String] Reference panel VCF with genotype information, for generation of .UD .Mean .Bed files[Optional]
--BamFile   [String] Bam or Cram file for the sample[Required]
--Reference [String] reference file[Required]
--Seed      [INT] Random number seed(default:12345)
--numPC     [INT] Number of Principal Components used in estimation
--fixPC     [String] Specify known PC coordinates for the sample, format(x.xxx|x.xxx)
--fixAlpha  [Float] Specify known contamination level
--asHeter   [Bool] Enable using hetergeneous model, in which intended sample and contaminating sample have different ancestries.(Recommended)
--knownAF   [String] A Bed file that provide known allele frequency for each marker, similar behaviour with VerifyBamID 1.0
```

## Generating your own resource files.

For producing customized resource files to be used as the input argument of verifyBamID2 generation, you need to start with a VCF file and FASTA formatted reference files. Please refer to the example below.
```
VerifyBamID --RefVCF ReferencePanel.vcf.gz --BamFile ./resource/test.bam --Reference ./resource/chr20.fa.gz

```

In the example above, the expected output file names will be ``ReferencePanel.vcf.gz.UD, ReferencePanel.vcf.gz.Mean, ReferencePanel.vcf.gz.Bed``

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
