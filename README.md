# radBARCODER

scripts to extract, align, and type mtDNA data from restriction site associated DNA sequenced on an [Illumina Machine](https://en.wikipedia.org/wiki/Illumina,_Inc.) with mitochondrial reference genomes of non-model species on linux/unix computers

* `radBARCODER bam2GENO`

* `radBARCODER aliGENO`

* `radBARCODER mkMETAGENO`

* `radBARCODER fltrGENOSITES`

---

## Quick Start

#### 1. Clone repo and put tools in path

```bash
# assuming OS = Ubuntu
git clone https://github.com/cbirdlab/radBARCODER.git
cd radBARCODER
sudo cp radB* /usr/local/bin
sudo cp *R /usr/local/bin
sudo chmod 555 /usr/local/bin/radBARCODER
sudo chmod 555 /usr/local/bin/fltrGENOSITES.R
sudo chmod 555 /usr/local/bin/consensusSEQ.R
```

#### 2. Install dependencies

These are the packages used by `radBARCODER`. If you install the newest versions, it will work.  More detailed instructions for installation can be found farther down in this document.

* [`parallel`](https://www.gnu.org/software/parallel/), [`bedtools`](https://github.com/arq5x/bedtools2/releases), [`samtools`](https://www.htslib.org/), [`bcftools`](https://samtools.github.io/bcftools/bcftools.html), [`pagan2`](http://wasabiapp.org/software/pagan/), [`mafft`](https://mafft.cbrc.jp/alignment/software/), [`seaview`](http://doua.prabi.fr/software/seaview) 

* [`R`](https://www.r-project.org/)

  * [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html), [stringr](https://cran.r-project.org/web/packages/stringr/index.html)

#### 3. Prep your data

All of the following files should be in a single directory:

* *`*bam` files*

  * Your NGS data should be in the form of binary alignment maps ([BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map)) that have been created using your NGS data and the reference mitochondrial genome. The `bam` file should also be indexed (see [tabix](https://github.com/samtools/tabix)).  1 `bam` and 1 `bam.bai` per individual. The following naming format is assumed: `Population*UniqueID*bam`.  If you want guidance on making the `bam` files, then scroll down to the detailed instructions that begin after the quick start instructions.
  
* *Sample classification files*

  * `radBARCODER` allows you to classify your samples into two groupings, independent of population classification, that will be used to construct consensus metagenomes.  The purpose of these files is to identify sequences as the expected target taxon or the unexpected taxon that showed up in your NGS data.  The files should be text and contain one ID per line, with the ID matching the pattern used in the `bam` files.  For example: `PopID_IndividualID`. Even if you do not have two groups of samples, you still need one file with the names of all the samples and a second file that is empty.  Naming format of the file, itself, is not important.

* *Mitochondrial genomes*

  * One mitochondrial genome should be selected to be the reference [FASTA](https://en.wikipedia.org/wiki/FASTA_format) with 1 sequence, not 1 sequence per gene region. This should be the file used to make the [binary alignment maps](https://en.wikipedia.org/wiki/Binary_Alignment_Map) from your NGS data. The reference can be one of the mitochondrial genomes below.

  * [optional] The mitochondrial genomes that you want to compare your data to should be in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format with 1 sequence, not 1 sequence per gene region. There should be a common pattern in name of these files, such as `*_mtGenomes.fasta` to allow `radBARCODER` to work correctly.  1 genome per file. *If you do include these additional genomes, be sure that your reference shares the common naming pattern.*
  
* *locus-specific fasta files*

  * [optional] If you also want to include incomplete mitochondrial sequences in the alignment, you need a single fasta file with as many sequences as you wish.  These sequences do not need to be aligned.


  
#### 4. `radBARCODER bam2GENO`: Convert `bam` files to mitochondrial genome sequences

Use `radBARCODER bam2GENO` to convert each of the `bam` files to a mitochondrial genome sequence.  All intermediate and final files are saved to `./out_bam2GENO`

```bash
# Name of reference mtGenome that you mapped your NGS reads to
REF=reference.Pfalc.mtGenome.fasta  

# String pattern match to the bam files for every individual (do not include a leading wildcard). Use same pattern match conventions as you would with `ls`
bamPATTERN=.Pfalc.mtGenome-RG.bam    

# number of processors to use for parallel operations
THREADS=32   

radBARCODER bam2GENO $REF $bamPATTERN $THREADS
```
  
#### 5. `radBARCODER aliGENO`: Align the mitochondrial genomes

Use `radBARCODER aliGENO` to align the mitchondrial genome sequences with either `pagan2` (`LONGALIGNMENT=FALSE`) or `mafft` (`LONGALIGNMENT=TRUE`). It pretty easy to modify the alignment method in the `align` function (found in `radBarcoder_functions.bash`). All intermediate and final files are saved to `./out_aliGENO`

```bash
REF=reference.Pfalc.mtGenome.fasta
bamPATTERN=.Pfalc.mtGenome-RG.bam
THREADS=32

# toggle between two alignment methods. FALSE is "better" than TRUE if it works.
LONGALIGNMENT=FALSE

# positions to align, refering to reference genome
POSITIONS=1-18000

# Name of locus or loci, will be included in output file names
LOCUS="PfalcMitoGenome"

# Name appended to beginning of output file names
PREFIX=paganAlign_

# Pattern match for all mitochondrial genomes to be included in alignment (use wildcards as necessary).  If you only have the reference genome, then use this `mtGenPATTERN=$REF` instead of the example pattern that follows
mtGenPATTERN="A*Genome.fasta"

# Name of fasta file containing partial mitochondrial sequences from GenBank
GENBANKFASTA=""

radBARCODER aliGENO $REF $bamPATTERN $THREADS $PREFIX $LOCUS $POSITIONS "$mtGenPATTERN" $LONGALIGNMENT $GENBANKFASTA
```

#### 6. `radBARCODER mkMETAGENO`: Make meta mitochondrial genomes

Use `radBARCODER mkMETAGENO` to create a consensus meta mitochondrial genomes for each population as well as two predefined groups of individuals in your NGS data. These are useful when you recover small portions of the mitochondrial genome in each individual. All intermediate and final files are saved to `./out_mkMETAGENO`

```bash
PREFIX=paganAlign_
LOCUS="PfalcMitoGenome"
THREADS=32

# list of sample names from individuals that are divergent from most of your samples to make a meta genome for 
nontargetIDs=$(cat Pproctozystron.txt)
nontargetNAME=Ppr

# list of sample names for the majority of your samples that genetically group together to make a meta genome for
targetIDs=$(cat Pfalcifer.txt)
targetNAME=Pfa

# a list of codes used to label population identity in the names of the `bam` files.  A meta genome will be made for each population from the individuals in the list of "targetIDs"
POPS=$(echo -e AtMk"\t"At"\t"Pk"\t"Kr"\t"St)

# for each position in the genome alignment, the minimum number of individuals having data that are required to include the consensus nucleotide call
cvgForCall=1

radBARCODER mkMETAGENO "$nontargetIDs" "$targetIDs" "$POPS" $PREFIX $LOCUS $THREADS $cvgForCall $nontargetNAME $targetNAME
```

#### 7. `radBARCODER fltrGENOSITES`: Selectively filter your final genome and meta genome alignments

Use `radBARCODER fltrGENOSITES` to filter genomes with more missing/ambiguous/indel base calls than specified with `PCT` then remove sites with missing/ambiguous/indel base calls.  Higher values of `PCT` retain more individuals and lower values retain more nuclotides. The result is a `fasta` alignment with only genomes and sites with A, C, T, or G nucleotide calls.  Output files include the `fasta` alignment, a `pdf` with descriptive plots, and a `csv` describing the reference genome positions retained.  Output files are saved to same directory as the input file (`FASTA`).

```bash
# name of file with mitochondrial genome or meta genome alignments
FASTA=paganAlign_ALL_masked_aligned_clean_PfalcMitoGenome.fasta

# keep nucleotides that occur in this percent of genomes and metagenomes in the alignment
PCT=99

radBARCODER fltrGENOSITES $FASTA $PCT
```


#### 8. Evolutionary reconstruction and haplotype networks

Use your favorite software to analyze and visualize the resulting alignments. I tend to use `popart` and `raxml`.


---

################################################################################

# DETAILED GUIDANCE

If you hit roadblocks in the quick start, more details are given below.  

We start by assuming that you have untrimmed `fq.gz` NGS files in a project directory, no `bam` files, and go from there.  I do assume that the `fq.gz` files have been demultiplexed but have not been trimmed or filtered.

All bash code assumes that your OS is [Ubuntu](https://ubuntu.com/).  


#### 1. Detailed Installation and Dependencies

We will use `dDocentHPC` to trim and map your files to the reference genome prior to using `radBARCODER` to generate mitochondrial genome alignments and metagenomes. Here, I assume that you will clone fresh copies of the `radBARCODER` and `dDocentHPC` repos into your project directory (so we avoid the vagaries of file permissions and the `$PATH`) and run the scripts directly rather than putting them into your `$PATH`.  

Clone the radBARCODER and dDocentHPC repos to your project dir:

```bash
# move to your directory for this project. replace "ProjectDir" with the path to the directory for this project
cd ProjectDir  

# clone repos to your project dir as follows
git clone https://github.com/cbirdlab/radBARCODER.git   
git clone https://github.com/cbirdlab/dDocentHPC.git
 
# move a few files to the project dir
cp dDocentHPC/config.4.all .
mkdir mkBAM
cp radBARCODER/radB* mkBAM
cp radBARCODER/*R mkBAM
```

The `config.4.all` file has all the settings for `dDocentHPC` and it is good practice to make a copy of this file in each dir from which it is run to serve as a record of what settings were used.

`radBARCODER` has fewer options and so does not include a config file. I do recommend that you save a bash script describing the settings and commands you ran, but I do not include that in the instructions below.

Assumed directory structure:

```
$ tree ../ProjectDir
ProjectDir
 ├──dDocentHPC
 ├──config.4.all
 ├──mkBAM
 │   ├──consensusSEQ.R
 │   ├──fltrGENOSITES.R
 │   ├──radBARCODER
 │   └──radBARCODER_functions.bash
 ├──pop1_ind1.F.fq.gz
 ├──pop1_ind1.R.fq.gz
 ...
 └──radBARCODER
```

Goto [`dDocentHPC`](https://github.com/cbirdlab/dDocentHPC) and find instructions to install all of the required software dependencies and clone the dDocentHPC repository. You will be asked to install [`dDocent`](https://www.ddocent.com, the program from which `dDocentHPC` was forked. There is a script that automatically installs the required software on your unix-based system. `dDocent` shares many similarities with `dDocentHPC` but the instructions here assume you are using `dDocentHPC`. You can run `dDocentHPC` on a workstation or HPC. 

If processing ddRAD libraries that are based on Peterson et al. (2012), I recommend adding 2 adapter sequences to the `trimmomatic` adapters file by overwriting the file `TruSeq3-PE-2.fa` that comes with `trimmomatic` with the modified `TruSeq3-PE-2.fa` file in the `radBARCODER` dir. You may also modify the `TruSeq3-PE-2.fa` file as necessary for your flavor of library prep.

```
# this will work on a workstation. on an HPC, run `dDocentHPC trimFQ` (see below) and view the output to see the path to the adapters file. 
# If trimmomatic has not been installed in /usr/local/bin, then the path below should be changed accordingly
cd ProjectDir
sudo cp radBARCODER/TruSeq3-PE-2.fa /usr/local/bin/adapters
```

`radBARCODER` has a few additional dependencies. Unfortunately, there is no installation script for them, but it is not difficult.  I provide some commands below which should work but it is up to you to find and update the URLs to the latest versions and make sure that the unzipped tarball dir names match the provided code.

* [`pagan2`](http://wasabiapp.org/software/pagan/) 

* [`mafft`](https://mafft.cbrc.jp/alignment/software/) 

* [`seaview`](http://doua.prabi.fr/software/seaview) 

* [`R`](https://www.r-project.org/)

  * [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html), [stringr](https://cran.r-project.org/web/packages/stringr/index.html)

```
# goto your downloads directory and update your apps (assuming Ubuntu or Debian OS)
cd ~/Downloads
sudo apt update
sudo apt upgrade

# install seaview (assuming Ubuntu or Debian OS)
wget http://doua.prabi.fr/software/seaview_data/seaview5-64.tgz
tar -xzf seaview5-64.tgz
# this could overwrite newer installations of clustalo, muscle. an alternative to this command is to add the pagan dir to the `$PATH`
sudo cp seaview/* /usr/local/bin

# install pagan2 (assuming Ubuntu or Debian OS)
wget http://wasabiapp.org/download/pagan/pagan2.linux64.20190829.tgz
tar -xzf pagan2.linux64.20190829.tgz
# this could overwrite newer installations of mafft, raxml. an alternative to this command is to add the pagan dir to the `$PATH`
sudo cp pagan/bin/* /usr/local/bin

# install mafft (assuming Ubuntu or Debian OS)
# note that pagan2 includes an old version of mafft and doing this could overwrite it.  I have not had any problems yet, but an alternative would be to add the pagan dir to the `$PATH`
wget https://mafft.cbrc.jp/alignment/software/mafft_7.471-1_amd64.deb
sudo dpkg -i mafft_7.471-1_amd64.deb

# See if you have R (ctrl-d to exit)
R

# if you do not have R, install R (assuming Ubuntu or Debian OS)
sudo apt -y install r-base

# alternatively, you could install the latest version of R from the GUI or following these CLI instructions https://cran.revolutionanalytics.com/doc/manuals/r-release/R-admin.html#Installing-R-under-Unix_002dalikes

# run R and install required packages
R
install.packages(c("seqinr", "stringr"))

# ctrl-d to exit R
```

---

#### 2. Detailed: Preparing your [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files & a reference mitchondrial genome 

Follow these steps to make mtGenomes from each individual in your RAD data set.  We use the [dDocentHPC](https://github.com/cbirdlab/dDocentHPC) pipeline for processing RAD data in unix-based computers.  It is assumed that your FASTQ files are minimally processed (demultiplexed with no quality trimming) gzipped and have the following naming convention : 

```
# files must end with F.fq.gz and R.fq.gz
# only 1 underscore should occur, and it should delimit the population idenity and the individual identity.  
# every individual must have a different identity
ProjectDir/Population_UniqueIndividualID.F.fq.gz
ProjectDir/Population_UniqueIndividualID.R.fq.gz
```

It is also assumed that you have a fully assembled mitochondrial genome saved as a [FASTA](https://en.wikipedia.org/wiki/FASTA) file.  The file should contain 1 sequence. You should name the reference mtGenome used for mapping sequence reads as follows:

```
# no more and no less than 3 periods should be used in the name and the * should be replaced with descriptive characters. I recommend species name and and genbank accession number.
ProjectDir/mkBAM/reference.*.*.fasta
```

Assumed directory structure:

```
$ tree ../ProjectDir
ProjectDir
 ├──dDocentHPC
 ├──config.4.all
 ├──mkBAM
 │   ├──consensusSeq.R
 │   ├──cullSeqs.R
 │   ├──maximizeBP.R
 │   ├──radBARCODER
 │   ├──radBARCODER_functions.bash
 │   └──reference.*.*.fasta
 ├──pop1_ind1.F.fq.gz
 ├──pop1_ind1.R.fq.gz
 ...
 └──radBARCODER
```

---

#### 3. `dDocentHPC trimFQmap`: Trim `fastq` Files for Mapping

Before running `dDocentHPC`, you should adjust the settings in the config file `config.4.all` as necessary.  `trimmomatic` is used to complete trimming which will remove low quality base calls, adapters, and reads that are too short after the removal of nucleotides.  The settings that affect reads trimmed for mapping to the mtGenome are labeled in the config as `mkBAM` signifying that these reads will be used to make the BAM files. Also, the `radBARCODER` repo has a modified adapter file which adds a ddRAD oligo that can make its way into the beginning of your sequence reads.  See installation section above.

```bash
nano config.4.all
```

relevant portion of `config.4.all` assuming ddRAD, 5bp barcodes, and 151 bp PE Illumina sequencing (settings that are `mkREF only` will not be used):

```
32              Number of Processors (Auto, 1, 2, 3, ..., n threads) cbirdq=40 normal=20
120G    Maximum Memory (1G,2G,..., 256G)  G=gigabytes
----------trimFQ: Settings for Trimming FASTQ Files---------------------------------------------------------------
146		trimmomatic MINLEN (integer, mkREF only)						Drop the read if it is below a specified length. Set to the length of the Read1 reads.
75		trimmomatic MINLEN (integer, mkBAM only)						Drop the read if it is below a specified length. Set to the minimum frag length you want mapped to the reference.
20		trimmomatic LEADING:<quality> (integer, mkBAM only)				Specifies the minimum quality required to keep a base.
15		trimmomatic TRAILING:<quality> (integer, mkREF only)			Specifies the minimum quality required to keep a base.
20		trimmomatic TRAILING:<quality> (integer, mkBAM only)			Specifies the minimum quality required to keep a base.
2		trimmomatic ILLUMINACLIP:<seed mismatches> (integer)			specifies the maximum mismatch count which will still allow a full match to be performed
30		trimmomatic ILLUMINACLIP:<palindrome clip thresh> (integer)		specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment
10		trimmomatic ILLUMINACLIP:<simple clip thresh> (integer)			specifies how accurate the match between any adapter etc. sequence must be against a read.
20		trimmomatic SLIDINGWINDOW:<windowSize> (integer)				specifies the number of bases to average across
20		trimmomatic SLIDINGWINDOW:<windowQuality> (integer)				specifies the average quality required.
0   trimmomatic CROP:<bp to keep> (integer, mkBAM only)    Trim read sequences down to this length. Enter 0 for no cropping
0		trimmomatic HEADCROP:<length> (integer, only Read1 for ezRAD)	The number of bases to remove from the start of the read. 0 for ddRAD, 5 for ezRAD
no		FixStacks (yes,no)   											Demultiplexing with stacks introduces anomolies.  This removes them.  
------------------------------------------------------------------------------------------------------------------
```

If your data is not ddRAD with 5bp barcodes or is not 151 bp PE, you should determine the lengths of reads in both your F and R `fq.gz` files and adjust `MINLEN` as you see fit.

Run dDocentHPC to trim the reads as follows:

```bash
# this could take a while and should probably be done on a powerful workstation (ex. >=16 threads, >=64 gb RAM)
# if you run out of memory, reduce the number of Processors specified in config.4.all and rerun
bash dDocentHPC/dDocentHPC.bash trimFQmap config.4.all
```

This will create a dir called `mkBAM` if it does not exist and it is populated with the trimmed `*R[12].fq.gz` files.


#### 4. `dDocentHPC mkBAM` and `dDocentHPC fltrBAM`: Map NGS data to mtDNA Genome then Filter


* obtain reference genome from [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)

  * reference genome should be a `fasta` formatted file and can be composed of 1, several, or all loci in the mtGenome
  
  * as an example, you could name reference genome as follows: `reference.GenusSpecies.GenBankAccession.fasta` 
  
* set cutoff in the `config.4.all` file to *_GenusSpecies_*

* set cutoff2 im the `config.4.all` file to *_GenBankAccession_*

Here is an example of the relevant portion of `config.4.all`:

```bash
----------mkREF: Settings for de novo assembly of the reference genome--------------------------------------------
PE              Type of reads for assembly (PE, SE, OL, RPE)                                    PE=ddRAD & ezRAD pairedend, non-overlapping reads;$
0.9             cdhit Clustering_Similarity_Pct (0-1)                                                   Use cdhit to cluster and collapse uniq rea$
Hspil           Cutoff1 (integer)                                                                                               Use unique reads t$
NC_023222               Cutoff2 (integer)                                                                                               Use unique$
0.05    rainbow merge -r <percentile> (decimal 0-1)                                             Percentile-based minimum number of seqs to assembl$
0.95    rainbow merge -R <percentile> (decimal 0-1)                                             Percentile-based maximum number of seqs to assembl$
------------------------------------------------------------------------------------------------------------------
```

Move to the mkBAM dir, map the reads to the mtGenome, and filter the resulting BAM files

```bash
# note that the paths assume the directory structure specified above
cd mkBAM

# this could take a while and should probably be done on a powerful workstation (ex. 16 threads, 64 gb RAM)
bash ../dDocentHPC/dDocentHPC.bash mkBAM ../config.4.all
bash ../dDocentHPC/dDocentHPC.bash fltrBAM ../config.4.all
```

This will create mildly filtered `*RG.bam` files for each individual. These alignment maps are used in downstream processing.  You should view the alignment maps with [IGV](https://software.broadinstitute.org/software/igv/download) or an equivalent `bam` viewer to ensure that mapping and filtering were successful.  Artifacts to look for are reads with too many SNPs (inappropriate alignment score threshold), large insertions at the beginning or ends of reads (adapter and barcode seqs not successfully trimmed), general sloppiness, etc.  If your reads were not originally 151 bp, you will probably need to change the alignment settings in the `config.4.all` file and rerun mapping.


#### 5. `radBARCODER bam2GENO` Create consensus sequences for each individual's reads mapped to the reference genome and mask areas with no coverage

*Dependencies*: [`parallel`](https://www.gnu.org/software/parallel/) [`bedtools`](https://github.com/arq5x/bedtools2/releases) [`samtools`](https://www.htslib.org/)  [`bcftools`](https://samtools.github.io/bcftools/bcftools.html) (fyi, all are required by `ddocent`, so you should have these installed if you made it to this step)

From here forward, you'll be running the `radBARCODER` scripts.  If you have not already, install required missing dependencies and clone the `radBARCODER` repo into your ProjectDir, and move the `radBARCODER*` and `*R` scripts to the mkBAM dir (see instructions above).

As a reminder, incontrast to the quick start section, here I describe running the `radBARCODER` scripts in your project directory. If you want to run them from the `$PATH`, you will have to modify the paths in the code blocks in the detailed instruction sections below. This is the expected dir structure after completing the previous step (some files and dirs created by trimming and mapping are omitted):

```
$ tree ../ProjectDir
ProjectDir
 ├──dDocentHPC
 ├──config.4.all
 ├──mkBAM
 │   ├──consensusSEQ.R
 │   ├──fltrGENOSITES.R
 │   ├──pop1_ind1.R1.fq.gz
 │   ├──pop1_ind1.R2.fq.gz
 ...
 │   ├──pop1_ind1.*.*.-RAW.bam
 │   ├──pop1_ind1.*.*.-RAW.bam.bai
 │   ├──pop1_ind1.*.*.-RG.bam
 │   ├──pop1_ind1.*.*.-RG.bam.bai
 ...
 │   ├──radBARCODER
 │   ├──radBARCODER_functions.bash
 │   └──reference.*.*.fasta
 ├──pop1_ind1.F.fq.gz
 ├──pop1_ind1.R.fq.gz
 ...
 └──radBARCODER
```

Update the following variable assignments and run `radBARCODER`:

```bash
# move to mkBAM dir.  replace "ProjectDir" as required given you dir structure.
cd ProjectDir/mkBAM

#Name of reference mtGenome
REF=reference.Pfalc.mtGenome.fasta  

#Pattern to id the bam files for each individual, this must be formatted as follows .CUTOFF1.CUTOFF2-RG.bam, where the CUTOFFs come from the dDocent config settings in step 2 above and are used in the name of the reference mtGenome
bamPATTERN=.Pfalc.mtGenome-RG.bam    

#number of processors to use for parallel operations
THREADS=8    

bash radBARCODER bam2GENO $REF $bamPATTERN $THREADS
```

Successful output looks like this:

```
#########################################################################
Sat 05 Sep 2020 12:51:38 AM CDT RUNNING radBARCODER BAM2GEN...
#########################################################################

Sat 05 Sep 2020 12:51:38 AM CDT VARIABLES READ IN:

the function that will be run FUNKTION=......bam2gen
the reference genome used to map the reads REF=...........reference.Pfalc.mtGenome.fasta
the ls pattern shared by all bam files bamPATTERN=.....Pfalc.mtGenome-RG.bam
the number of cpu cores for the task THREADS=.......4
the characters added to every file created PREFIX=........
the name of the locus or loci LOCUS=.........
the nucleotide positions in the reference genome to consider POSITIONS=.....
the ls pattern shared by all mtGenomes that will be aligned mtGenPATTERN=..
the aligner that will be used LONGALIGNMENT=.
the GenBank sequences that should also be aligned GENBANK=.......


Sat 05 Sep 2020 12:51:38 AM CDT 190 SAMPLES BEING PROCESSED:
AtMk_Pfa052 AtMk_Pfa055 AtMk_Pfa077 AtMk_Pfa086 AtMk_Pfa089 AtMk_Pfa095 At_Pfa034 At_Pfa035 At_Pfa036 At_Pfa037 At_Pfa038 At_Pfa039 At_Pfa040 At_Pfa041 At_Pfa043 At_Pfa044 At_Pfa045 At_Pfa046 At_Pfa048 At_Pfa049 At_Pfa050 At_Pfa054 At_Pfa056 At_Pfa058 At_Pfa059 At_Pfa060 At_Pfa061 At_Pfa062 At_Pfa063 At_Pfa065 At_Pfa072 At_Pfa074 At_Pfa075 At_Pfa076 At_Pfa078 At_Pfa079 At_Pfa080 At_Pfa081 At_Pfa082 At_Pfa083 At_Pfa084 At_Pfa085 At_Pfa087 At_Pfa088 At_Pfa091 At_Pfa092 At_Pfa094 Kr_Pfa001 Kr_Pfa003 Kr_Pfa004 Kr_Pfa005 Kr_Pfa007 Kr_Pfa008 Kr_Pfa010 Kr_Pfa011 Kr_Pfa012 Kr_Pfa013 Kr_Pfa014 Kr_Pfa015 Kr_Pfa016 Kr_Pfa018 Kr_Pfa020 Kr_Pfa021 Kr_Pfa022 Kr_Pfa024 Kr_Pfa025 Kr_Pfa027 Kr_Pfa029 Kr_Pfa030 Kr_Pfa031 Kr_Pfa032 Kr_Pfa033 Kr_Pfa035 Kr_Pfa037 Kr_Pfa039 Kr_Pfa040 Kr_Pfa041 Kr_Pfa042 Kr_Pfa043 Kr_Pfa044 Kr_Pfa045 Kr_Pfa046 Kr_Pfa047 Kr_Ppr002 Kr_Ppr006 Kr_Ppr009 Kr_Ppr017 Kr_Ppr019 Kr_Ppr023 Kr_Ppr026 Kr_Ppr028 Kr_Ppr034 Kr_Ppr036 Kr_Ppr038 Pk_Pfa321 Pk_Pfa322 Pk_Pfa324 Pk_Pfa325 Pk_Pfa326 Pk_Pfa327 Pk_Pfa328 Pk_Pfa329 Pk_Pfa330 Pk_Pfa331 Pk_Pfa332 Pk_Pfa333 Pk_Pfa334 Pk_Pfa335 Pk_Pfa336 Pk_Pfa337 Pk_Pfa338 Pk_Pfa339 Pk_Pfa340 Pk_Pfa341 Pk_Pfa343 Pk_Pfa346 Pk_Pfa347 Pk_Pfa348 Pk_Pfa349 Pk_Pfa351 Pk_Pfa352 Pk_Pfa353 Pk_Pfa356 Pk_Pfa357 Pk_Pfa358 Pk_Pfa359 Pk_Pfa360 Pk_Pfa361 Pk_Pfa362 Pk_Pfa363 Pk_Pfa364 Pk_Pfa365 Pk_Pfa366 Pk_Pfa367 Pk_Pfa370 Pk_Pfa371 Pk_Pfa372 Pk_Pfa373 Pk_Pfa375 Pk_Pfa376 Pk_Pfa377 Pk_Pfa380 St_Pfa001 St_Pfa002 St_Pfa003 St_Pfa004 St_Pfa005 St_Pfa006 St_Pfa007 St_Pfa008 St_Pfa009 St_Pfa010 St_Pfa011 St_Pfa012 St_Pfa013 St_Pfa014 St_Pfa015 St_Pfa016 St_Pfa018 St_Pfa019 St_Pfa021 St_Pfa022 St_Pfa023 St_Pfa024 St_Pfa025 St_Pfa026 St_Pfa027 St_Pfa028 St_Pfa029 St_Pfa030 St_Pfa031 St_Pfa032 St_Pfa033 St_Pfa034 St_Pfa035 St_Pfa036 St_Pfa037 St_Pfa040 St_Pfa041 St_Pfa042 St_Pfa043 St_Pfa044 St_Pfa045 St_Pfa046 St_Pfa047 St_Pfa048 St_Pfa049 St_Pfa050 St_Pfa051 St_Ppr017

[mpileup] 1 samples in 1 input files
[mpileup] maximum number of reads per input file set to -d 30000
Lines   total/split/realigned/skipped:  904/0/0/0
Applied 37 variants
[mpileup] 1 samples in 1 input files
[mpileup] maximum number of reads per input file set to -d 30000
Lines   total/split/realigned/skipped:  734/0/0/0
Applied 51 variants
...
[mpileup] 1 samples in 1 input files
[mpileup] maximum number of reads per input file set to -d 30000
Lines   total/split/realigned/skipped:  1087/0/0/0
Applied 31 variants

#########################################################################
Sat 05 Sep 2020 12:51:46 AM CDT radBARCODER BAM2GEN COMPLETED
#########################################################################
```

This should result in a `vcf.gz` and a `masked_consensus.fasta` for every individual.

Hard coded stringencies are:

* depth of coverage >= 1
  * modify by changing `'\t0$'` in `radBARCODER_functions.bash`
    * example:  `'\t9$'` will mask all positions with < 10 reads

* base call phred quality >= 20
  * modify `-Q 20` in line beginning with `bcftools mpileup` in `radBARCODER_functions.bash`

* mapping quality >= 30
  * modify `-q 30` in line beginning with `bcftools mpileup` in `radBARCODER_functions.bash`
  
 Note that heterozygous positions are set to default to the reference allele. This behavior can be modified in `radBARCODER_functions.bash` at the line beginning with `bcftools consensus`. 

Updated directory structure:

```
$ tree ../../ProjectDir
ProjectDir
 ├──dDocentHPC
 ├──config.4.all
 ├──mkBAM
 │   ├──consensusSEQ.R
 │   ├──fltrGENOSITES.R
 │   ├──out_bam2GENO
 ...
 │   │  ├──pop1_ind1.*.*-RG_masked_*vcf.gz
 │   │  ├──pop1_ind1.*.*-RG_masked_consensus.fasta
 ...
 │   │  └──all.*.*-RG_maksed_consensus.fasta
 │   ├──pop1_ind1.R1.fq.gz
 │   ├──pop1_ind1.R2.fq.gz
 ...
 │   ├──pop1_ind1.*.*-RAW.bam
 │   ├──pop1_ind1.*.*-RAW.bam.bai
 │   ├──pop1_ind1.*.*-RG.bam
 │   ├──pop1_ind1.*.*-RG.bam.bai
 ...
 │   ├──radBARCODER
 │   ├──radBARCODER_functions.bash
 │   └──reference.*.*.fasta
 ├──pop1_ind1.F.fq.gz
 ├──pop1_ind1.R.fq.gz
 ...
 └──radBARCODER
```


#### 6. `radBARCODER aliGENO`: Select All or a Portion of the Genomes and Align Them Among Individuals

*Dependencies*: [`pagan`](http://wasabiapp.org/software/pagan/) [`mafft`](https://mafft.cbrc.jp/alignment/software/) [`seaview`](http://doua.prabi.fr/software/seaview) 

If you have additional sequences from GenBank that you would like to include with this alignment, you should download them and save into your 'mkBAM' dir. Additional mtGenome sequences should be saves as FASTA files and renamed to have a common format of your choosing. A FASTA file for targeted locus sequences can also be downloaded and included in the alignment.  See the specification of user-defined variables `mtGenPATTERN` and `GENBANKFASTA` below. 

Note that seaview is only used to convert from `fasta` to `nexus` format, so if you don't have it installed, you can manually convert the `fasta` to `nexus`. Also, while the `pagan2` precompiled `tar.gz` does have `mafft`, it is not complete and you should install the complete `mafft` if you are aligning very long sequences, otherwise you will get an error when setting `LONGALIGNMENT=TRUE`

You can specify which positions to target (for specific loci and genes) to make alignments, including disjunct positions. For example, if you want to specify positions 1-10, then:

```bash
POSITIONS=1-10
```

If you want to align only positions 1-4 and 6-10, then:

```bash
POSITIONS=1-4,6-10
```

Refer to the mtGenomes annotation for the positions of particular loci of interest.

The other issue is which aligner to use.  I have tried `clustalw`, `clustalo`, `mafft`, and `pagan2`.  I've found `pagan2` to be superior in that it almost never needs to be aligned by eye to clean up mistakes.  If it does, you should scrutinize your `bam` files with `IGV` and possibly rerun previous steps with new settings.  The default behavior is to run `pagan2` for alignment.  However, if you run into problems, potentially due to sequences being very long, then you can try `mafft` with options set for very long sequences as follows:

```bash
LONGALIGNMENT=TRUE   #use mafft instead of pagan2
```

Update the following variable assignments and run `radBARCODER`:

```bash
REF=reference.Pfalc.mtGenome.fasta
bamPATTERN=.Pfalc.mtGenome-RG.bam
THREADS=4

# positions from the reference mtGenome and BAM files to include in the alignment
POSITIONS=1-18000

# the name of the locus or loci specified in POSITIONS
LOCUS="PfalcMitoGenome"

#prefix on output files 
PREFIX=paganAlign_

# a pattern using wildcards that matches all mtGenome files other than the reference mtGenome
mtGenPATTERN="A*Genome.fasta"

# the name of a fasta file that has sequences to include in the alignment
GENBANKFASTA=""

# specify whether pagan (FALSE) or mafft (TRUE) is used to align.  pagan is better
LONGALIGNMENT=FALSE

radBARCODER aliGENO $REF $bamPATTERN $THREADS $PREFIX $LOCUS $POSITIONS "$mtGenPATTERN" $LONGALIGNMENT $GENBANKFASTA
```

Successful output looks like this:

```
#########################################################################
Sat 05 Sep 2020 01:09:40 AM CDT RUNNING radBARCODER ALIGN...
#########################################################################

Sat 05 Sep 2020 01:09:40 AM CDT VARIABLES READ IN:

the function that will be run FUNKTION=......align
the reference genome used to map the reads REF=...........reference.Pfalc.mtGenome.fasta
the ls pattern shared by all bam files bamPATTERN=.....Pfalc.mtGenome-RG.bam
the number of cpu cores for the task THREADS=.......4
the characters added to every file created PREFIX=........paganAlign_
the name of the locus or loci LOCUS=.........PfalcMitoGenome
the nucleotide positions in the reference genome to consider POSITIONS=.....1-18000
the ls pattern shared by all mtGenomes that will be aligned mtGenPATTERN=..A*Genome.fasta
the aligner that will be used LONGALIGNMENT=.FALSE
the GenBank sequences that should also be aligned GENBANK=.......


Sat 05 Sep 2020 01:09:40 AM CDT 190 SAMPLES BEING PROCESSED:
AtMk_Pfa052 AtMk_Pfa055 AtMk_Pfa077 AtMk_Pfa086 AtMk_Pfa089 AtMk_Pfa095 At_Pfa034 At_Pfa035 At_Pfa036 At_Pfa037 At_Pfa038 At_Pfa039 At_Pfa040 At_Pfa041 At_Pfa043 At_Pfa044 At_Pfa045 At_Pfa046 At_Pfa048 At_Pfa049 At_Pfa050 At_Pfa054 At_Pfa056 At_Pfa058 At_Pfa059 At_Pfa060 At_Pfa061 At_Pfa062 At_Pfa063 At_Pfa065 At_Pfa072 At_Pfa074 At_Pfa075 At_Pfa076 At_Pfa078 At_Pfa079 At_Pfa080 At_Pfa081 At_Pfa082 At_Pfa083 At_Pfa084 At_Pfa085 At_Pfa087 At_Pfa088 At_Pfa091 At_Pfa092 At_Pfa094 Kr_Pfa001 Kr_Pfa003 Kr_Pfa004 Kr_Pfa005 Kr_Pfa007 Kr_Pfa008 Kr_Pfa010 Kr_Pfa011 Kr_Pfa012 Kr_Pfa013 Kr_Pfa014 Kr_Pfa015 Kr_Pfa016 Kr_Pfa018 Kr_Pfa020 Kr_Pfa021 Kr_Pfa022 Kr_Pfa024 Kr_Pfa025 Kr_Pfa027 Kr_Pfa029 Kr_Pfa030 Kr_Pfa031 Kr_Pfa032 Kr_Pfa033 Kr_Pfa035 Kr_Pfa037 Kr_Pfa039 Kr_Pfa040 Kr_Pfa041 Kr_Pfa042 Kr_Pfa043 Kr_Pfa044 Kr_Pfa045 Kr_Pfa046 Kr_Pfa047 Kr_Ppr002 Kr_Ppr006 Kr_Ppr009 Kr_Ppr017 Kr_Ppr019 Kr_Ppr023 Kr_Ppr026 Kr_Ppr028 Kr_Ppr034 Kr_Ppr036 Kr_Ppr038 Pk_Pfa321 Pk_Pfa322 Pk_Pfa324 Pk_Pfa325 Pk_Pfa326 Pk_Pfa327 Pk_Pfa328 Pk_Pfa329 Pk_Pfa330 Pk_Pfa331 Pk_Pfa332 Pk_Pfa333 Pk_Pfa334 Pk_Pfa335 Pk_Pfa336 Pk_Pfa337 Pk_Pfa338 Pk_Pfa339 Pk_Pfa340 Pk_Pfa341 Pk_Pfa343 Pk_Pfa346 Pk_Pfa347 Pk_Pfa348 Pk_Pfa349 Pk_Pfa351 Pk_Pfa352 Pk_Pfa353 Pk_Pfa356 Pk_Pfa357 Pk_Pfa358 Pk_Pfa359 Pk_Pfa360 Pk_Pfa361 Pk_Pfa362 Pk_Pfa363 Pk_Pfa364 Pk_Pfa365 Pk_Pfa366 Pk_Pfa367 Pk_Pfa370 Pk_Pfa371 Pk_Pfa372 Pk_Pfa373 Pk_Pfa375 Pk_Pfa376 Pk_Pfa377 Pk_Pfa380 St_Pfa001 St_Pfa002 St_Pfa003 St_Pfa004 St_Pfa005 St_Pfa006 St_Pfa007 St_Pfa008 St_Pfa009 St_Pfa010 St_Pfa011 St_Pfa012 St_Pfa013 St_Pfa014 St_Pfa015 St_Pfa016 St_Pfa018 St_Pfa019 St_Pfa021 St_Pfa022 St_Pfa023 St_Pfa024 St_Pfa025 St_Pfa026 St_Pfa027 St_Pfa028 St_Pfa029 St_Pfa030 St_Pfa031 St_Pfa032 St_Pfa033 St_Pfa034 St_Pfa035 St_Pfa036 St_Pfa037 St_Pfa040 St_Pfa041 St_Pfa042 St_Pfa043 St_Pfa044 St_Pfa045 St_Pfa046 St_Pfa047 St_Pfa048 St_Pfa049 St_Pfa050 St_Pfa051 St_Ppr017


Sat 05 Sep 2020 01:09:40 AM CDT EXTRACTING POSITIONS 1-18000 FOR ALIGNMENT...

Sat 05 Sep 2020 01:09:40 AM CDT REMOVING INDIVIDUALS WITH NO NUCLEOTIDES CALLED...

Sat 05 Sep 2020 01:09:40 AM CDT GATHERING ALL mtGENOMEs WITH PATTERN=A*Genome.fasta AND EXTRACTING POSITIONS 1-18000

Sat 05 Sep 2020 01:09:41 AM CDT ADDING reference.Pfalc.mtGenome.fasta SEQUENCES FROM GENBANK...

Sat 05 Sep 2020 01:09:41 AM CDT CONCATENATING FASTAs...

Sat 05 Sep 2020 01:09:41 AM CDT ALIGNING SEQUENCES WITH pagan2...

Sat 05 Sep 2020 01:09:41 AM CDT ALIGNING RAD DATA TO MITOGENOMES...

PAGAN2 v.1.53 (29 August, 2019). (C) 2010-2019 by Ari Löytynoja <ari.loytynoja@gmail.com>.
 This program is provided "as-is", with NO WARRANTY whatsoever; this is a development version
 and may contain bugs.

The analysis started: Sat Sep  5 01:09:41 2020
Alignment file: ref.fas
Guidetree file: ref.tre

The analysis finished: Sat Sep  5 01:09:42 2020
Total time used by PAGAN: 0.796133 wall sec, 0.459117 cpu sec.


Sat 05 Sep 2020 01:09:52 AM CDT Making sure all individuals have same number of nucleotides plus indels
     There are 16593 nucleotides in the longest sequence
     953 indels being added to short sequence

#########################################################################
Sat 05 Sep 2020 01:10:49 AM CDT radBARCODER ALIGN COMPLETED
#########################################################################
```

The final genome alignments are in the files named `${PREFIX}_ALL_masked_aligned_clean_${LOCUS}.` in the `out_aliGENO` dir, and will be the newest files in the directory (use `ls -ltrh out_aliGENO`) to view. Many intermediate files are also found in the same dir.  `${PREFIX}` and `${LOCUS}` represent the values save to those variables in the code above.

Updated directory structure:

```
$ tree ../../ProjectDir
ProjectDir
 ├──dDocentHPC
 ├──config.4.all
 ├──mkBAM
 │   ├──consensusSEQ.R
 │   ├──fltrGENOSITES.R
 │   ├──out_aliGENO
 ...
 │   │  ├──${PREFIX}_ALL_masked_aligned_clean_${LOCUS}.fasta
 │   │  └──${PREFIX}_ALL_masked_aligned_clean_${LOCUS}.nex
 │   ├──out_bam2GENO
 ...
 │   │  ├──pop1_ind1.*.*-RG_masked_*vcf.gz
 │   │  ├──pop1_ind1.*.*-RG_masked_consensus.fasta
 ...
 │   │  └──all.*.*-RG_maksed_consensus.fasta
 │   ├──pop1_ind1.R1.fq.gz
 │   ├──pop1_ind1.R2.fq.gz
 ...
 │   ├──pop1_ind1.*.*-RAW.bam
 │   ├──pop1_ind1.*.*-RAW.bam.bai
 │   ├──pop1_ind1.*.*-RG.bam
 │   ├──pop1_ind1.*.*-RG.bam.bai
 ...
 │   ├──radBARCODER
 │   ├──radBARCODER_functions.bash
 │   └──reference.*.*.fasta
 ├──pop1_ind1.F.fq.gz
 ├──pop1_ind1.R.fq.gz
 ...
 └──radBARCODER
```

It is important to check the alignment by eye and edit as necessary or adjust upstream settings. I recommend [`seaview`](http://doua.prabi.fr/software/seaview) for this, but any alignment viewer will work. In the example data set, which has a lot of missing data, I did not have adjust the alignment at all, but mileage may vary.


#### 7. `radBARCODER mkMETAGEN` Make Meta-Genomes

*Dependencies*: `R` (`seqinr`, `stringr`) 

If you didn't have much luck extracting the same portions of the mtGenome among individuals in steps 1-6, you can make meta mitochondrial genomes from groups of individuals and align those using `mkMETAGEN` 

This function will make consensus sequences for each sample category following the dDocent naming convention (`PopulationID_IndividualID`), but you need to specify the the population ids as described below.  The genesis of radBARCODER was trying to figure out what an unexpected population partition was, so it is also assumed that a subset of individuals will be identified at "nonTarget". A text file with one id (`PopulationID_IndividualID`) per line can be used for this as shown below.  The remaining individuals belonging to the majority or targeted taxon should be listed similarly in a separate file. A metagenome will be created for the individuals in each of these files. The two files with the sample names must not contain any of the same individuals.  If there are no "nontarget" individuals, then replace the line below with `nontargetIDs=""`.  

Abbreviated contents of `Pproctozystron.txt`, which contains nontarget taxon that was not expected:

```
St_Ppr017
Kr_Ppr002
Kr_Ppr006
Kr_Ppr009
Kr_Ppr017
```

Abbreviated contents of `Pfalcifer.txt`, which contains nontarget taxon that was not expected:

```
St_Pfa047
St_Pfa048
St_Pfa049
St_Pfa050
St_Pfa051
```

You can also create a meta genome for each population, as long as it is coded into the names of the files. In the example code below, the populations are `At_`, `AtMk_`, `Pk_`, `Kr_`, `St_`. Note the that `_` will prevent individuals in population `AtMk` from being included in `At`.  The `\t` are tab delimiters. Spaces in place of the tabs will probably work also, but no commas should be used. These population identifiers are used to select the individals from the targeted taxon listed in `$targetIDs`. 

`cvgForCall` will determine the minimum read depth required to make a consensus base call. Base calling is performed by `consensusSeq.R` if you want to modify.

```bash
PREFIX=paganAlign_
LOCUS="PfalcMitoGenome"
THREADS=32

# list of sample names from individuals that are divergent from most of your samples to make a meta genome for 
nontargetIDs=$(cat Pproctozystron.txt)
nontargetNAME=Ppr

# list of sample names for the majority of your samples that genetically group together to make a meta genome for
targetIDs=$(cat Pfalcifer.txt)
targetNAME=Pfa

# a list of codes used to label population identity in the names of the `bam` files.  A meta genome will be made for each population from the individuals in the list of "targetIDs"
POPS=$(echo -e AtMk"\t"At"\t"Pk"\t"Kr"\t"St)

# for each position in the genome alignment, the minimum number of individuals having data that are required to include the consensus nucleotide call
cvgForCall=1

bash radBARCODER mkMETAGENO "$nontargetIDs" "$targetIDs" "$POPS" $PREFIX $LOCUS $THREADS $cvgForCall $nontargetNAME $targetNAME
```

See below for expected output metagenome files. Several other intermediate files will also be save to `out_metaGENO`. `${PREFIX}` and `${LOCUS}` represent the values save to those variables in the code above.

Updated directory structure:

```
$ tree ../../ProjectDir
ProjectDir
 ├──dDocentHPC
 ├──config.4.all
 ├──mkBAM
 │   ├──consensusSEQ.R
 │   ├──fltrGENOSITES.R
 │   ├──out_aliGENO
 ...
 │   │  ├──${PREFIX}_ALL_masked_aligned_clean_${LOCUS}.fasta
 │   │  └──${PREFIX}_ALL_masked_aligned_clean_${LOCUS}.nex
 │   ├──out_bam2GENO
 ...
 │   │  ├──pop1_ind1.*.*-RG_masked_*vcf.gz
 │   │  ├──pop1_ind1.*.*-RG_masked_consensus.fasta
 ...
 │   │  └──all.*.*-RG_maksed_consensus.fasta
 │   ├──out_metaGENO
 ...
 │   │  ├──${PREFIX}_ALL_masked_aligned_clean_metageno_${LOCUS}.fasta
 │   │  ├──${PREFIX}_ALL_masked_aligned_clean_metageno_${LOCUS}.nex
 │   │  ├──${PREFIX}_NonTargetTaxon_masked_aligned_metageno_${LOCUS}.fasta
 │   │  ├──${PREFIX}_NonTargetTaxon_masked_aligned_metageno_${LOCUS}.nex
 │   │  ├──${PREFIX}_pop1-POP_masked_aligned_clean_metageno_${LOCUS}.fasta
 │   │  ├──${PREFIX}_pop1-POP_masked_aligned_clean_metageno_${LOCUS}.nex
...
 │   │  ├──${PREFIX}_POPS_masked_aligned_clean_metageno_${LOCUS}.fasta
 │   │  ├──${PREFIX}_POPS_masked_aligned_clean_metageno_${LOCUS}.nex
 │   │  ├──${PREFIX}_TargetTaxon_masked_aligned_metageno_${LOCUS}.fasta
 │   │  ├──${PREFIX}_TargetTaxon_masked_aligned_metageno_${LOCUS}.nex
 │   │  ├──${PREFIX}_TAXA_masked_aligned_clean_metageno_${LOCUS}.fasta
 │   │  └──${PREFIX}_TAXA_masked_aligned_clean_metageno_${LOCUS}.nex
 │   ├──pop1_ind1.R1.fq.gz
 │   ├──pop1_ind1.R2.fq.gz
 ...
 │   ├──pop1_ind1.*.*-RAW.bam
 │   ├──pop1_ind1.*.*-RAW.bam.bai
 │   ├──pop1_ind1.*.*-RG.bam
 │   ├──pop1_ind1.*.*-RG.bam.bai
 ...
 │   ├──radBARCODER
 │   ├──radBARCODER_functions.bash
 │   └──reference.*.*.fasta
 ├──pop1_ind1.F.fq.gz
 ├──pop1_ind1.R.fq.gz
 ...
 └──radBARCODER
```


#### 8. `radBARCODER fltrGENOSITES`  Selectively filter your final genome and meta genome alignments

*Dependencies*: `R` (`seqinr`, `stringr`) 

This function will call `fltrGENOSITES.R` which is included in the `radBARCODER` repo.  Make sure it is in your `mkBAM` directory

Filter genomes with more missing/ambiguous/indel base calls than specified with `PCT` then remove sites with missing/ambiguous/indel base calls.  Higher values of `PCT` retain more individuals and lower values retain more nuclotides. The result is a `fasta` alignment with only genomes and sites with A, C, T, or G nucleotide calls.  Output files include the `fasta` alignment, a `pdf` with descriptive plots, and a `csv` describing the reference genome positions retained.  Output files are saved to same directory as the input file (`FASTA`).

Update the following variable assignments and run `radBARCODER`:

```bash
# path to aligned (meta) genome file
FASTA=out_aliGENO/paganAlign_ALL_masked_aligned_clean_PfalcMitoGenome.fasta

# remove (meta) genomes with > `$PCT` missing/ambiguous/indel nucleotide calls.
PCT=99

radBARCODER fltrGENOSITES $FASTA $PCT

PCT=95
radBARCODER fltrGENOSITES $FASTA $PCT
PCT=75
radBARCODER fltrGENOSITES $FASTA $PCT
PCT=50
radBARCODER fltrGENOSITES $FASTA $PCT
PCT=25
radBARCODER fltrGENOSITES $FASTA $PCT
PCT=5
radBARCODER fltrGENOSITES $FASTA $PCT
PCT=1
radBARCODER fltrGENOSITES $FASTA $PCT
```

Example successful output from 1 run looks like this (while some errors show up here, bute these are ok because they are not related to the `cullGENO` function or are normal output that I have not suppressed (yet):

```
#########################################################################
Sat 05 Sep 2020 01:19:12 AM CDT RUNNING radBARCODER fltrGENOSITES...
#########################################################################

Sat 05 Sep 2020 01:19:12 AM CDT VARIABLES READ IN:

the function that will be run FUNKTION=......fltrGENOSITES
the reference genome used to map the reads REF=...........paganAlign_ALL_masked_aligned_clean_PfalcMitoGenome.fasta
the ls pattern shared by all bam files bamPATTERN=....99
the number of cpu cores for the task THREADS=.......
the characters added to every file created PREFIX=........
the name of the locus or loci LOCUS=.........
the nucleotide positions in the reference genome to consider POSITIONS=.....
the ls pattern shared by all mtGenomes that will be aligned mtGenPATTERN=..
the aligner that will be used LONGALIGNMENT=.
the GenBank sequences that should also be aligned GENBANK=.......

Sat 05 Sep 2020 01:19:12 AM CDT 1 SAMPLES BEING PROCESSED:

Sat 05 Sep 2020 01:19:12 AM CDT Be sure to check the fasta alignment by eye as small errors do occur
[1] "paganAlign_ALL_masked_aligned_clean_PfalcMitoGenome"
[2] "99"
null device
          1

#########################################################################
Sat 05 Sep 2020 01:19:13 AM CDT radBARCODER fltrGENOSITES COMPLETED
#########################################################################
```

You can also check for proper output files to confirm success:

```
$ ls -tr | tail -n3
paganAlign_ALL_masked_aligned_clean_PfalcMitoGenome_99.fasta
missingCalls99.pdf
paganAlign_ALL_masked_aligned_clean_PfalcMitoGenome_99.nex
```

You might see an error message referring to `max` and `min` if `PCT` is too large for your data, but the output should still be created.

The output files include a `pdf` with informative figures describing the (meta) genomes, a `csv` with the positions retained by the filter, and `fasta/nexus` files with the aligned (meta) genomes.  You can cross reference the positions in the `csv` against the genome annotation to determine the loci retained after filtering.



#### 9. Make network with `PopArt` 

[`PopArt`](https://github.com/jessicawleigh/popart-current), or your favorite network program, can now be used to create a network from the file.  `PopArt` automatically removes positions and sequences with poor coverage, so it's very convenient to apply to the file at this point.  [Precompiled, but outdated versions of PopArt](http://popart.otago.ac.nz/index.shtml)

The `nex` (meta) genome alignments can be read directly into `popart`.


#### 10. Determine the best evolutionary model with [`modeltest-ng`](https://github.com/ddarriba/modeltest)

The `fasta` (meta) genome alignments can be read directly into `modeltest-ng`. If you are reading this, then I figure you might need some help installing `modeltest-ng` because I recall having dependency issues that were not handled in their documentation. 

if you install `modeltest-ng`, make sure you have the following dependencies:

```bash
sudo apt-get install flex bison libgmp3-dev
```


#### 11. Reconstruct evolutionary histories

[raxml](https://github.com/amkozlov/raxml-ng) can be used to resconstruct evolutionary histories from the `fasta` formatted (meta) genome alignments.  There are also [web servers](https://raxml-ng.vital-it.ch/#/) that work well for small data sets that do not require large amounts of bootstraps.

if you install `raxml`, make sure you have the following dependencies:

```bash
sudo apt-get install flex bison libgmp3-dev
```

