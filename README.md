# VerifyBamID

Motivation: Detecting and estimating inter-sample DNA contamination became a crucial quality assessment step to ensure high quality sequence reads and reliable downstream analysis. Across many existing models, allele frequency usually is used to calculate prior genotype probability. Lack of accurate assignment of allele frequency could result in underestimation of contamination level.
Results: We applied our method to 1000 Genomes datasets by simulating contamination levels from 1% to 20% and comparing the contamination estimates obtain from different methods. When using pooled allele frequencies, as opposed to population-specific allele frequencies, we observed that the contamination levels are underestimated by 20%, 40%, 51%, and 73% for CEU, YRI, FIN, and CHS populations, respectively. Using our new method, the underestimation bias was reduced to 2-5%.


## Installation

  - mkdir build
  - cd build
  - cmake ..
  - make
  - make test

## Usage
```
VerifyBamID --UDPath ./resource/hapmap_3.3.b37.dat.UD --BamFile ./resource/test.bam --BedPath ./resource/choose.bed --MeanPath ./resource/hapmap_3.3.b37.dat.mu --Reference ./resource/chr20.fa.gz
```
```
--UDPath    [String] UD matrix file from SVD result of genotype matrix[Required]
--MeanPath  [String] Mean matrix file of genotype matrix[Required]
--BamFile   [String] Bam or Cram file for the sample[Required]
--BedPath   [String] Bed file for markers used in this analysis,(chr\tpos-1\tpos\trefAllele\taltAllele)[Required]
--Reference [String] reference file[Required]
--Seed      [INT] Random number seed(default:12345)
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
