# tHapMix - lightweight and fast simulation of tumour whole-genome sequencing data

## Table of Contents
 - [Citation](#citation)
 - [Description](#description)
 - [License](#license)
 - [System Requirements](#system-requirements)
 - [Installation and dependencies](#installation-and-dependencies)
  - [Installing from binaries](#installing-from-binaries)
  - [Installing from source code](#installing-from-source-code)
 - [Phased Haplotype BAMs](#phased-haplotype-bams)
 - [Usage](#usage)
 - [Examples](#examples)
 - [References](#references)

## Citation
tHapMix: simulating tumour samples through haplotype mixtures 
Sergii Ivakhno; Camilla Colombo; Stephen Tanner; Philip Tedder; Stefano Berri; Anthony J. Cox
Bioinformatics 2016; doi: 10.1093/bioinformatics/btw589.

## Description
tHapMix is a lightweight simulation engine that generates somatic copy number
and SNV variants on different cancer evolutionary trees to produce diverse
polyclonal tumour samples. tHapMix utilizes phased sequencing data, where
individual germline variants are split into haplotypes providing
chromosome-specific coverage and allele tracks. By sampling from such high
coverage tracks one can simulate copy number variants, SNVs and other features of
tumour samples such as ploidy, purity and heterogeneity. In particular, the following  features can be altered to recreate specific types of tumour evolution:

1. Number of clones 

2. Clonal evolutionary tree structure

3. Percentage of major clone 

4. Percentage of private (heterogeneous) variants

tHapMix can be used for software benchmarking, model training and in-silico analysis of tumour evolutionary mechanisms.

Examples of tHapMix data can be downloaded from https://basespace.illumina.com/projects/30822795 .
 
## License

tHapMix source code is provided under the [GPLv3 license] (LICENSE.txt).


##System Requirements: 

tHapMix uses the pyFlow workflow engine for parallel execution. It requires at least 8G of RAM for the whole-genome simulation, however for the high-coverage whole-genome simulation a multi-core / cluster system configuration is advisable. We tested and report runtime for the following configurations (80x whole-genome BAM):
- 4.5h on a Linux node with 8 CPUs and 60G RAM
- 70 minutes on a Linux node with 32 CPUs and 120G RAM


## Installation and dependencies
For all modes of installation, *bedtools* and *samtools* packages should be installed and must be present in *PATH* environment, see http://bedtools.readthedocs.org/en/latest/content/installation.html and http://www.htslib.org/ for details. tHapMix has been tested with *Python* versions *2.6* and *2.7* and should also run on *Python 2.8*. tHapMix was not tested with *Python 3*. Example installation for Ubuntu will look like:
- sudo apt-get install zlib1g-dev # for samtools 
- sudo apt-get install samtools
- sudo apt-get install bedtools

### Installing from binary archives
The easiest way to install tHapMix is by using provided binary distribution. To install tHapMix into ```$THAPMIXDIR```, just download a tar.gz archive for a corresponding version from releases github page (e.g. tHapMix-1.0.tar.gz) and run 

```bash
pip install -r requirements.txt
pip install HapMix-*.tar.gz -t $THAPMIXDIR --install-option="--install-data=$THAPMIXDIR"
```

This will also install all the dependencies (python-dev might need to be installed as well). After this tHapMix can be launched by a simple command

```bash
python $THAPMIXDIR/bin/submitHapMix.py -h
```

### Installing from source code

Below we list two installation options: via virtualenv and Docker. Docker installation also includes a small *demo example*. 

##### Installing with virtualenv

* Install openssl prerequisite for Python
```bash
wget https://www.openssl.org/source/openssl-1.0.1j.tar.gz
tar xzf openssl-1.0.1j.tar.gz
cd openssl-1.0.1j
CFLAG=-FPIC
./config -fPIC shared zlib-dynamic
make
make install
cd ..
```
* Install Python
```bash
wget https://www.python.org/ftp/python/2.7.5/Python-2.7.5.tgz
tar xzf Python-2.7.5.tgz
cd Python-2.7.5
./configure --disable-ipv6 --with-pydebug
make
make install
cd ..
```
* Install and build virtualenv
```bash
wget http://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz --no-check-certificate
tar xzf virtualenv-1.11.6.tar.gz
cd virtualenv-1.11.6
python virtualenv.py DESTDIR
```
* Install Python packages
```bash
DESTDIR/bin/pip2.7 install cython
DESTDIR/bin/pip2.7 install numpy
DESTDIR/bin/pip2.7 install pysam
DESTDIR/bin/pip2.7 install bx-python
```
* Use Python executable to run tHapMix
```
DESTDIR/bin/python2.7
```
##### Installing with Docker / Docker test run
The Dockerfile which is included in the top level directory of the tHapMix source distribution can be used to build an image as and run test example as:
```bash
0. Prerequisites
     - Example requires a minimum of 8G RAM and 40G of disk space 
     - For AmazonAWS, using east-coast instances (US East, N. Virginia) will enable faster file transfer
     - Examples below were tested on Ubuntu 
     - sudo apt-get install git -y
     - sudo apt-get install docker.io -y

1. Installation
     - git clone https://github.com/Illumina/tHapMix
     - cd HapMix ## this example presumes /tmp/tHapMix working directory 
     - sudo docker build   -t thapmix ./ # build docker image 

2. Usage 
     - Register (if needed) for BaseSpace and open https://basespace.illumina.com/s/Tty7T2ppH3Tr
     
     - The required files are available in this subfolder https://basespace.illumina.com/analyses/30880851/files/28484471?projectId=18065049
     
     - copy sorted_haplotypeC_chr21.bam, sorted_haplotypeC_chr22.bam, sorted_haplotypeD_chr21.bam, sorted_haplotypeD_chr22.bam plus corresponding bam indices into  /tmp/tHapMix/haplotyped_bams
     
         * You can use BaseMount (https://basemount.basespace.illumina.com) as a command line interface for the file transfer
         * If you chose to use BaseMount, the files will be in: <mount point>"/Projects/HiSeq 2000: TruSeq PCR-Free (Platinum Genomes)/AppResults/HaplotypeBamsNA12882/Files"
         * If you run, mkdir /tmp/Base_Space; basemount /tmp/Base_Space, <mount point> will equal /tmp/Base_Space
         
     - sudo docker run -v /tmp/tHapMix/haplotyped_bams:/opt/thapmix/haplotyped_bams -i -t thapmix /bin/bash # launch container
     
     - cd bin # run HapMix test example withing docker, see Examples section for more details
     
     - ./submitHapMix.py -c ../example/example.json -t ../example/example.bed -o ../haplotyped_bams/out (will take ~10min on a node with 2 cores)
     
     - exit
     
     - sorted_example_som_var.bam bam will be available in /tmp/tHapMix/haplotyped_bams/out/
     
     - The bam or the whole folder could be transferred to permanent folder.
```
where /tmp is the path to haplotype bam files. This will create a bash session from which tHapMix can be run.

## Phased Haplotype BAMs
Phased Haplotype SNV and indel variants from Platinum Genomes (PG) project (http://www.illumina.com/platinumgenomes/) comprising an extended 17-member family have been used to split aligned reads into separate haplotypes for each chromosome. The phased haplotypes were created by using MERLIN linkage software [1] on informative SNV calls generated by GATK [2]. PG sample NA12882 sequenced to a 200x depth on a HiSeq 2000 system was used a source of aligned reads. The variant selection went through a number of stringent criteria including: 

1)     Variant calls have no Mendelian inconsistencies within a pedigree.

2)     Variants are called in all replicates.

The splitting procedure generated a separate bam file for each of chromosomal haplotypes (i.e. 48 files taking a total of 300G of disk space). These can be downloaded from BaseSpace: https://basespace.illumina.com/s/Tty7T2ppH3Tr (exact folder location: https://basespace.illumina.com/analyses/30880851/files/28484471?projectId=18065049). A **BaseMount** command line interface for BaseSpace could be used for file downloading: details accessible from here **https://basemount.basespace.illumina.com**. 

It is possible to regenerate your own haplotype bams by using `split_by_haplotype.py` under `scripts/haplobams/` directory. To use it install the following Python package dependencies with `pip install [packagename]: pysam, numpy, HTSeq, subprocess` and make sure HTSlib is available in the `$PATH` (follow installation instructions in http://www.htslib.org/download). Once installed, `split_by_haplotype.py` will require the following files:
- phased SNV calls for Platinum Genome samples: available for download from tHapMix release 
- phased indel calls for Platinum Genome samples: available for download from tHapMix release 
- chromosome-level fasta genome reference: i.e., from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/
- NA12882 and NA12878 200x bam files: raw fastq files for the sample are available for download from European Nucleotide Archive under the accession number ERP001775.

To create haplotype bams run the following command for each chromosome bam 

```split_by_haplotype.py [-h] -s SAMPLE_NAME -c CHR_NO -b BAM_FILE -r REF_FILE -n SNP_FILE -i INDEL_FILE  [-x CHR_PREFIX]```.

For example, executing 

```python split_by_haplotype.py -s NA12878 -c 20 -b NA12878_chr20_small.bam -r chr20.fa -n IlluminaPlatinumSNPs-PG-staging-V4-all.vcf.gz -i IlluminaPlatinumINDELs-PG-staging-V4-all.vcf.gz```

will produce haplotype bams `haplotypeC_chr20.bam` and `haplotypeD_chr20.bam`.


## Usage
Run
```bash
python $THAPMIXDIR/bin/submitHapMix.py -h
```
for full usage instructions.

tHapMix uses the following input files 

1.  Either ground truth bed file  with known copy numbers per haplotype or
parameters for simulation
 * config/.bed files provide locations of CNVs while config/.vcf file - SNVs 
 * scripts/getCNVFeatures.py can be used to generate distribution of CNV features (size, copy numbers and genotypes) for use in simulation from a set of previously called CNVs in a VCF format
 * config generated by getCNVFeatures.py can found in config/varsimconfig/config_tumor_truth_file.json
 * to simulate CNV de-novo run  ```submitHapMix.py -b config/varsimconfig/config_tumor_truth_file.json -c config/config_fixed_single_level.json -o OUTDIR```
 * varsimconfig/config_tumor_truth_file.json can be re-generated from a set of VCF files by running ```scripts/getCNVFeatures.py -o config/varsimconfig/ -t 4```
 * to fix random number generation ```--seed``` option should be used

2.  Haplotype bams. 
 * haplotyped_bam_dir in in json configuration file: ("/path/to/directory/with/bams", see "Phased Haplotype BAMs" section above for details on downloading required bams)

3.  Simulates cancer sample properties in json configuration file (example
    values in parentheses, file can be found in config/.json)

 * purity: purity of the sample (60.0)
 * num_clones: total number of clones in a sample (4)
 * perc_majority_clone: percentage abundance of the major clone, abundance of 
   the rest is evenly split (50.0)
 * perc_heterogeneous_variants: percentage of variants that are not shared by
   all clones (100.0)
 * tree_format: "random_binary", "random_single_level" or custom (["A",["B","C"], "E"])
 * majority_clone: root clone ("A")
 * ploidy_depth: (80) - for diploid genome

Other configuration parameters for config/.json include:
 * bam_depth: (174)
 * mutate_sv: (true)
 * somatic_vcf: ("/path/to/somatic/snvs/file.vcf",) - the currently supported SNV VCF format is the one generated by Strelka somatic SNV caller, https://sites.google.com/site/strelkasomaticvariantcaller/. Example is shown in config/somatic.snvs.vcf
 * hg_file: ("/path/to/genomes/human.hg19.genome")
 

##Output: 
After tHapMix run output_dir contains the following files and subdirectories
 * Simulated tumour bam file
 * _clone{A}.bed - location and copy numbers of variants provate to clone {A}, etc.
 * sim_param.json - file with simulation parameters
 * somatic.snvs_clone{A}.vcf - SNV variants provate to clone {A}, etc.
 * pyflow.data/logs/ - log files
 * pyflow.data/state/ - commands used

##Examples
###Small example
Use config file in example subfolder for a simplified test run. Only chromosomes 21 and 22 will be simulated. For this you will need to create a haplotyped_bams directory directly under tHapMix installation and copy haplotype bams and bam indices for chromosome 21 and 22 into it. Then run
```
cd bin
submitHapMix.py -c ../example/example.json -t ../example/example.bed  -o ./examplerun
```
### Full demo
This demo will show how to perform full tHapMix simulation followed by CNV identification and analysis of the results. Supporting files for this demo can also be found in example subfolder.


##### Haplotype bams
NA12878 and NA12882 haplotype bams can be downloaded from BaseSpace: https://basespace.illumina.com/s/Tty7T2ppH3Tr (exact folder location: https://basespace.illumina.com/analyses/30880851/files/28484471?projectId=18065049). A **BaseMount** command line interface for BaseSpace could be used for file downloading (via e.g. cp): installation and usage instructions are available at **https://basemount.basespace.illumina.com**. Briefly, to mount ```$basespace``` folder type ```basemount $basespace```.

##### tHapMix simulation
For this simulation we will use an annotated truth set for pseudo-haploid genome (```GeL004flt.bed```) that has a ploidy of 1.6 and approximately 200 copy number changes and no LOH. NA12882 haplotype bams will be used as read input.  In addition, we will simulate the following sample and genome properties as encoded in sim_param.json file (entire configuration file for the demo can be found in ):
* Percentage of heterogeneous variants set to 70% ("var_het": 70.0)
* Purity set to 80% ("purity": 80.0)
* Sample with two clones ("num_clones": 2,)
* Abundance of major clone 60% ("perc_maj": 50.0)
With these input parameters the tHapMix simulation can be launched as 

```bash
python bin/submitHapMix.py -t config/GeL004flt.bed -c config/simGeL004flt_n2_m60_v70.json -o outdir -m local
```

##### Canvas CNV calling
tHapMix run will generate bam file in the output directory specified. We could assess simulated copy number changes by using a CNV calling tool. We use Canvas in this demo. Full documentation can be accessed at https://github.com/Illumina/Canvas. Here we show example command line as applied to tHapMix output. An already tHapMix pre-generated example bams are available for download from https://basespace.illumina.com/s/zkWcr91RIdLG: Analyses/bams subdirectory contains simulated tumour and normal bams, while Analyses/utils - SNV calls in vcf formal in the normal sample.
```
Canvas.exe Somatic-WGS 
--bam = path to tHapMix bam
--sample-name = simGeL004flt_n2_m70_v70
--reference = UCSC/hg19/Annotation/Canvas/kmer.fa
--genome-folder = UCSC/hg19/Sequence/WholeGenomeFasta
--filter-bed = UCSC/hg19/Annotation/Canvas/filter13.bed' (flagged region to be skipped, i.e. centromeres)
--output 
--b-allele-vcf normal.vcf.gz (vcf containing SNV b-allele sitesonly sites with PASS in the filter column will be used)
```


##References
[1] Abecasis,G. et al. (2002) Merlin - rapid analysis of dense genetic maps using sparse gene flow trees. Nat Genet., 20, 97-101.

[2] McKenna,A. et al. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res., 20, 1297-1303.
