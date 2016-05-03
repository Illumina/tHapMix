# HapMix - lightweight and fast simulation of tumour whole-genome sequencing data
## Description
HapMix is a lightweight simulation engine that generates somatic copy number
and SNV variants on different cancer evolutionary trees to produce diverse
polyclonal tumour samples. HapMix utilizes phased sequencing data, where
individual germline variants are split into haplotypes providing
chromosome-specific coverage and allele tracks. By sampling from such high
coverage tracks one can simulate copy number variants, SNVs and other features of
tumour samples such as ploidy, purity and heterogeneity. In particular, the following  features can be altered to recreate specific types of tumor evolution:

1. Number of clones 

2. Clonal evolutionary tree stracture

3. Percentage of major clone 

4. Percentage of private (heterogeneous) variants

HapMix can be used for software benchmarking, model training and in-silico analysis of tumour evolutionary mechanisms.
 
## License

HapMix source code is provided under the [GPLv3 license] (LICENSE.txt).

## Installation and dependencies
Bedtools package should be installed and must be present in *PATH* environment: this could be accomplished with package managers http://bedtools.readthedocs.org/en/latest/content/installation.html . 

HapMix relies on pyflow workflow engine (https://github.com/Illumina/pyflow) for tasks management and integrates samtools for bam file manipulation. To install these modules run *make* in HapMix *redist* subdirectory.

In addition pysam and bx-python needs to be available in PYTHOPATH, for installation instructions visit https://bitbucket.org/james_taylor/bx-python/wiki/HowToInstall and https://github.com/pysam-developers/pysam 

Below we list two installation options: via virtualenv and Docker. Docker installation also includes a small *demo example*. 

### Installing with virtualenv
Tested Python versions for HapMix so far include 2.7.5. To ease dependency maintenance a recommended way of installing HapMix is via virtualenv. The following steps with install all necessary packages. 

* Install openssl prerequisite for Python
```
wget https://www.openssl.org/source/openssl-1.0.1j.tar.gz
tar xzf openssl-1.0.1j.tar.gz
cd openssl-1.0.1j
CFLAG=-FPIC
./config -fPIC --prefix=/illumina/thirdparty/openssl --openssldir=/illumina/thirdparty/openssl shared zlib-dynamic
make
make install
cd ..
```
* Install Python
```
wget https://www.python.org/ftp/python/2.7.5/Python-2.7.5.tgz
tar xzf Python-2.7.5.tgz
cd Python-2.7.5
./configure --prefix=/illumina/thirdparty/python/python-2.7.5-dev --disable-ipv6 --with-pydebug
make
make install
cd ..
```
* Install and build virtualenv
```
wget http://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz --no-check-certificate
tar xzf virtualenv-1.11.6.tar.gz
cd virtualenv-1.11.6
/illumina/thirdparty/python/python-2.7.5-dev/bin/python virtualenv.py /illumina/thirdparty/python/python-2.7.5-dev/virtualenv
```
* Install Python packages
```
/illumina/thirdparty/python/python-2.7.5-dev/virtualenv/bin/pip2.7 install cython
/illumina/thirdparty/python/python-2.7.5-dev/virtualenv/bin/pip2.7 install numpy
/illumina/thirdparty/python/python-2.7.5-dev/virtualenv/bin/pip2.7 install pysam
/illumina/thirdparty/python/python-2.7.5-dev/virtualenv/bin/pip2.7 install bx-python
```
* Use Python executable to run HapMix
```
/illumina/thirdparty/python/python-2.7.5-dev/virtualenv/bin/python2.7
```
### Installing with Docker / Docker test run
The Dockerfile which is included in the top level directory of the HapMix source distribution can be used to build an image as and run test example as:
```
0. Prerequisites
     - Example requires a minimum of 8G RAM and 40G of disk space 
     - For AmazonAWS, using east-coast instances (US East, N. Virginia) will enable faster file transfer
     - Examples below were tested on Ubuntu 
     - sudo apt-get install git -y
     - sudo apt-get install docker.io -y

1. Installation
     - git clone https://git.illumina.com/Bioinformatics/HapMix
     - cd HapMix ## this example presumes /tmp/HapMix working directory 
     - sudo docker build   -t hapmix ./ # build docker image 

2. Usage 
     - Register (if needed) for BaseSpace and open https://basespace.illumina.com/s/Tty7T2ppH3Tr
     
     - The required files are available in this subfolder https://basespace.illumina.com/analyses/30880851/files/28484471?projectId=18065049
     
     - copy sorted_haplotypeC_chr21.bam, sorted_haplotypeC_chr22.bam, sorted_haplotypeD_chr21.bam, sorted_haplotypeD_chr22.bam plus corresponding bam indices into  /tmp/HapMix/haplotyped_bams
     
         * You can use BaseMount (https://basemount.basespace.illumina.com) as a command line interface for the file transfer
         * If you chose to use BaseMount, the files will be in: <mount point>"/Projects/HiSeq 2000: TruSeq PCR-Free (Platinum Genomes)/AppResults/HaplotypeBamsNA12882/Files"
         * You you run, mkdir /tmp/Base_Space; basemount /tmp/Base_Space, <mount point> will equal /tmp/Base_Space
         
     - sudo docker run -v /tmp/HapMix/haplotyped_bams:/opt/hapmix/haplotyped_bams -i -t hapmix /bin/bash # launch container
     
     - cd bin # run HapMix test example withing docker, see Examples section for more details
     
     - ./submitHapMix.py -c ../example/example.json -t ../example/example.bed -o ../haplotyped_bams/out (will take ~30min on a node with 2 cores)
     
     - exit
     
     - sorted_example_som_var.bam bam will be available in /tmp/HapMix/haplotyped_bams/out/
     
     - The bam or the whole folder could be transferred to permanent folder.
```
where /tmp is the path to haplotype bam files. This will create bash session from which HapMix can be run.

### Phased Haplotype BAMs
Phased Haplotype SNV and indel variants from Platinum Genomes (PG) project (http://www.illumina.com/platinumgenomes/) comprising an extended 17-member family have been used to split aligned reads into separate haplotypes for each chromosome. The phased haplotypes were created by using MERLIN linkage software [1] on informative SNV calls generated by GATK [2]. PG sample NA12882 sequenced to a 200x depth on a HiSeq 2000 system was used a source of aligned reads. The variant selection went through a number of stringent criteria including: 

1)     Variant calls have no Mendelian inconsistencies within a pedigree.

2)     Variants are called in all replicates.

The splitting procedure generated a separate bam file for each of chromosomal haplotypes (i.e. 48 files taking a total of 300G of disk space). These can be downloaded from BaseSpace: https://basespace.illumina.com/s/Tty7T2ppH3Tr (exact folder location: https://basespace.illumina.com/analyses/30880851/files/28484471?projectId=18065049). A BaseMount command line interface for BaseSpace could be used for file downloading: details accessible from here https://basemount.basespace.illumina.com . 


## Usage
HapMix uses the following input files 

1.  Either ground truth file bed with known copy numbers per haplotype or
parameters for simulation

2.  Haplotype bams.

3.  Simulates cancer sample properties in json configuration file (example
    values in parentheses)

 * purity: purity of the sample (60.0)
 * num_clones: total number of clones in a sample (4)
 * perc_majority_clone: percentage abundance of the major clone, abundance of 
   the rest is evenly split (50.0)
 * perc_heterogeneous_variants: percentage of variants that are not shared by
   all clones (100.0)
 * tree_format: "random_binary", "random_single_level" or custom (["A",["B","C"], "E"])
 * majority_clone: root clone ("A")
 * ploidy_depth: (60)

Other configuration parameters include:
* haplotyped_bam_dir: ("/path/to/directory/with/bams", see "Phased Haplotype BAMs" section above for details on downloading required bams)
* bam_depth: (174)
* mutate_sv: (true)
* somatic_vcf: ("/path/to/somatic/snvs/file.vcf",)
* hg_filev: ("/path/to/genomes/human.hg19.genome")
* output_dir: ("/path/to/result/directory")
* output_bam_file: ("root_name_for_files")* 

##Output: 
After HapMix run output_dir contains the following files and subdirectories
* Simulated tumour bam file
* _clone{A}.bed - location and copy numbers of variants provate to clone {A}, etc.
* sim_param.json - file with simulation parameters
* somatic.snvs_clone{A}.vcf - SNV variants provate to clone {A}, etc.
* pyflow.data/logs/ - log files
* pyflow.data/state/ - commands used

##System Requirements: 
HapMix uses pyflow workflow engine for parallel execution. It requires at least 4G of RAM for the whole-genome simulation .

##Examples
Use config file in example subfolder for a simplified test run. Only chromosomes 21 and 22 will be simulated. For this you will need to create a haplotyped_bams directory directly under HapMix installation and copy haplotype bams and bam indices for chromosome 21 and 22 into it. Then run
```
cd bin
submitHapMix.py -c ../example/example.json -t ../example/example.bed  -o ./examplerun
```
##References
[1] Abecasis,G. et al. (2002) Merlin - rapid analysis of dense genetic maps using sparse gene flow trees. Nat Genet., 20, 97-101.

[2] McKenna,A. et al. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res., 20, 1297-1303.

