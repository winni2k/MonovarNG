## MonovarNG - Monovar but better
## Overview ##
**MonovarNG** is a single nucleotide variant (SNV) detection and genotyping algorithm for single-cell DNA sequencing data. It takes a list of bam files as input and outputs a vcf file containing the detected SNVs. This improved version of Monovar is written in C++ and achieves 100x speedup over the original Monovar.

## Dependencies ##
* [Boost C++](http://boost.org)
* [Htslib](http://htslib.org)
* [Cmake](http://cmake.org)

## Installation ##
Install dependencies
```
sudo apt-get install git g++ cmake libboost-all-dev zlib1g-dev libbz2-dev liblzma-dev
```
Install htslib
```
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar -xvf htslib-1.9.tar.bz2
cd htslib-1.9
./configure
make
sudo make install
```
Build Monovar
```
cd
git clone https://github.com/tianyi5309/MonovarNG
cd MonovarNG
cmake .
make
```
Add Monovar to path
```
export PATH=$PATH:$CURR_DIR/bin/monovar
```

## Usage ##
The program requires multiple bam files. The bam files should be sorted by coordinates. The raw sequence reads in .fastq format should be aligned to a reference genome with the help of an aligner program (e.g., BWA ([http://bio-bwa.sourceforge.net/]())). Aligner like BWA generates sam files containing aligned reads. The sam files can be converted to compressed bam files using ```samtools view``` command (see Samtools manual for details [http://www.htslib.org/doc/samtools.html]()). 


```
samtools mpileup -BQ0 -d10000 -f ref.fa -q 40 -b filenames.txt | monovar ref.fa filenames.txt compiled.pl output.vcf [-patm]
```
The arguments of Monovar are as follows:

```
-t: Threshold to be used for variant calling (Recommended value: 0.05)
-p: Offset for prior probability for false-positive error (Recommended value: 0.002)
-a: Offset for prior probability for allelic drop out (Default value: 0.2)
-m: Number of threads to use in multiprocessing (Default value: 1)
```
We recommend using cutoff 40 for mapping quality when using ```samtools mpileup```. To use the probabilistic realignment for the computation of Base Alignment Quality, drop the ```-B``` while running ```samtools mpileup```.
