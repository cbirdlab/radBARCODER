# radBARCODER

scripts to extract, align, and type mtDNA data from restriction site associated DNA sequenced on an [Illumina Machine](https://en.wikipedia.org/wiki/Illumina,_Inc.) with mitochondrial reference genomes of non-model species

---

## Installation and Dependencies

It is up to you how to handle the radBARCODER and dDocentHPC scripts, but here I assume that you will clone fresh copies of the two repos into your project directory and run the scripts directly rather than putting them into your `$PATH`.  Clone the radBARCODER and dDocentHPC repos to your project dir:

```bash
# move to your directory for this project. replace "ProjectDir" with the path to the directory for this project
cd ProjectDir  

# clone repos to your project dir as follows
git clone https://github.com/cbirdlab/radBARCODER.git   
git clone https://github.com/cbirdlab/dDocentHPC.git
 
# set up dir
cp dDocentHPC/config.4.all .
mkdir mkBAM
cp radBARCODER/*bash mkBAM
cp radBARCODER/*R mkBAM
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
 │   ├──radBARCODER.bash
 │   └──radBarcoder_functions.bash
 ├──pop1_ind1.F.fq.gz
 ├──pop1_ind1.R.fq.gz
 ...
 └──radBARCODER
```

Goto [dDocentHPC](https://github.com/cbirdlab/dDocentHPC) and find instructions to install all of the required software dependencies and clone the dDocentHPC repository. There is a script that automatically installs the software on your unix-based system. dDocentHPC was forked from [dDocent](https://www.ddocent.com) and shares many similarities but the instructions here assume you are using dDocentHPC. You can run dDocentHPC on a workstation or HPC. It is up to you whether you put the `dDocentHPC.bash` script into your `$PATH` or run it directly from the repo.  I usually clone a fresh copy to the top level of a project directory and execute it directly with `bash`.

If processing ddRAD libraries that are based on Peterson et al. (2012), I recommend adding 2 adapter sequences to the `trimmomatic` adapters file by overwriting the file `TruSeq3-PE-2.fa` which comes with `trimmomatic` with the modified `TruSeq3-PE-2.fa` file in the `radBARCODER` dir. You may also modify the `TruSeq3-PE-2.fa` file as necessary for your flavor of library prep.

```
# this will work on a workstation. on an HPC, run `dDocentHPC trimFQ` (see below) and view the output to see the path to the adapters file
cd ProjectDir
sudo cp radBARCODER/TruSeq3-PE-2.fa /usr/local/bin/adapters
```

`radBARCODER` has a few additional dependencies. Unfortunately, there is no installation script for them but it is not difficult.  I provide some commands below which should work but it is up to you to find and update the URLs to the latest versions and make sure that the unzipped tarball dir names match the provided code.

* [`pagan2`](http://wasabiapp.org/software/pagan/) 

* [`mafft`](https://mafft.cbrc.jp/alignment/software/) 

* [`seaview`](http://doua.prabi.fr/software/seaview) 

* [`R`](https://www.r-project.org/)

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
```

---

## Preparing your [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files & a reference mitchondrial genome 

Follow these steps to make mtGenomes from each individual in your RAD data set.  We use the [dDocentHPC](https://github.com/cbirdlab/dDocentHPC) pipeline for processing RAD data in unix-based computers.  It is assumed that your FASTQ files are minimally processed (demultiplexed with no quality trimming) gzipped and have the following naming convention : 

```
# files must end with [FR].fq.gz
# only 1 underscore should occur, and it should delimit the population idenity and the individual identity.  
# every individual must have a different identity
ProjectDir/Population_UniqueIndividualID.F.fq.gz
ProjectDir/Population_UniqueIndividualID.R.fq.gz
```

It is also assumed that you have a fully assembled mitochondrial genome saved as a [FASTA](https://en.wikipedia.org/wiki/FASTA) file. You should name the reference mtGenome used for mapping sequence reads as follows:

```
# no more and no less than 3 periods should be used in the name and the * should be replaced with descriptive characters.
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
 │   ├──radBARCODER.bash
 │   ├──radBarcoder_functions.bash
 │   └──reference.*.*.fasta
 ├──pop1_ind1.F.fq.gz
 ├──pop1_ind1.R.fq.gz
 ...
 └──radBARCODER
```

---

## Quick Start

#### 1. Trim `fastq` files for mapping: [dDocentHPC trimFQmap](https://github.com/cbirdlab/dDocentHPC)

Before running `dDocentHPC`, you should adjust the settings in the config file `config.4.all` as necessary.  `trimmomatic` is used to complete trimming which will remove low quality base calls, adapters, and reads that are too short after the removal of nucleotides.  The settings that affect reads trimmed for mapping to the mtGenome are labeled in the config as `mkBAM` signifying that these reads will be used to make the BAM files.

```bash
nano config.4.all
```

relevant portion of `config.4.all`:

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

Run dDocentHPC to trim the reads as follows:

```bash
# this could take a while and should probably be done on a powerful workstation (ex. >=16 threads, >=64 gb RAM)
# if you run out of memory, reduce the number of Processors specified in config.4.all and rerun
bash dDocentHPC/dDocentHPC.bash trimFQmap config.4.all
```

This will create a dir called `mkBAM` if it does not exist and it is populated with the trimmed `*R[12].fq.gz` files.


#### 2. Map `fastq` to mtDNA genome using [dDocentHPC mkBAM](https://github.com/cbirdlab/dDocentHPC)


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

This will create mildly filtered `RG.bam` files for each individual. These alignment maps are used in downstream processing.  You should view the alignment maps with [IGV](https://software.broadinstitute.org/software/igv/download) or an equivalent bam viewer to ensure that mapping and filtering were successful.  Artifacts to look for are reads with too many SNPs (inappropriate alignment score threshold), large insertions at the beginning or ends of reads (adapter and barcode seqs not successfully trimmed), etc.  If your reads were not originally 150 bp, you will probably need to change the alignment settings in the `config.4.all` file and rerun mapping.


#### 3. Create consensus sequences for each individual's reads mapped to the reference genome and mask areas with no coverage using `bam2fasta`

*Dependencies*: [`parallel`](https://www.gnu.org/software/parallel/) [`bedtools`](https://github.com/arq5x/bedtools2/releases) [`samtools`](https://www.htslib.org/)  [`bcftools`](https://samtools.github.io/bcftools/bcftools.html) (fyi, all are required by `ddocent`, so you should have these if you made it to this step)

From here forward, you'll be running the `radBARCODER.bash` scripts.  If you have not already, install required missing dependencies and clone the `radBARCODER` repo into your ProjectDir, and move the `radBARCODER` `*bash` and `*R` scripts to the mkBAM dir (see instructions above).

As a reminder, this is the expected dir structure after completing the previous step (some files and dirs created by trimming and mapping are omitted):

```
$ tree ../ProjectDir
ProjectDir
 ├──dDocentHPC
 ├──config.4.all
 ├──mkBAM
 │   ├──consensusSeq.R
 │   ├──cullSeqs.R
 │   ├──maximizeBP.R
 │   ├──pop1_ind1.R1.fq.gz
 │   ├──pop1_ind1.R2.fq.gz
 ...
 │   ├──pop1_ind1.*.*.-RAW.bam
 │   ├──pop1_ind1.*.*.-RAW.bam.bai
 │   ├──pop1_ind1.*.*.-RG.bam
 │   ├──pop1_ind1.*.*.-RG.bam.bai
 ...
 │   ├──radBARCODER.bash
 │   ├──radBarcoder_functions.bash
 │   └──reference.*.*.fasta
 ├──pop1_ind1.F.fq.gz
 ├──pop1_ind1.R.fq.gz
 ...
 └──radBARCODER
```

Update the following variable assignments and run `radBARCODER`:

```bash
#Name of reference mtGenome
REF=reference.Pfalc.mtGenome.fasta  

#Pattern to id the bam files for each individual, this must be formatted as follows .CUTOFF1.CUTOFF2-RG.bam, where the CUTOFFs come from the dDocent config settings in step 2 above and are used in the name of the reference mtGenome
bamPATTERN=.Pfalc.mtGenome-RG.bam    

#number of processors to use for parallel operations
THREADS=8    

bash radBARCODER.bash bam2fasta $REF $bamPATTERN $THREADS
```

Successful output looks like this:

```
#########################################################################
Sat 05 Sep 2020 12:51:38 AM CDT RUNNING radBARCODER BAM2FASTA...
#########################################################################

Sat 05 Sep 2020 12:51:38 AM CDT VARIABLES READ IN:

the function that will be run FUNKTION=......bam2fasta
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
Sat 05 Sep 2020 12:51:46 AM CDT radBARCODER BAM2FASTA COMPLETED
#########################################################################
```

This should result in a `vcf.gz` and a `masked_consensus.fasta` for every individual. Note that heterozygous positions are set to default to the reference allele. This behavior can be modified in `radBarcodder_functions.bash` at the line beginning with `bcftools consensus`. 

Hard coded strigencies are:

* depth of coverage >= 1

* base call phred quality >= 20

* mapping quality >= 30


#### 4. Select a portion of the genomes and `align` it across individuals

*Dependencies*: [`pagan`](http://wasabiapp.org/software/pagan/) [`mafft`](https://mafft.cbrc.jp/alignment/software/) [`seaview`](http://doua.prabi.fr/software/seaview) 

If you have additional sequences from GenBank that you would like to include with this alignment, you should download them and save into your 'mkBAM' dir. Additional mtGenome sequences should be saves as FASTA files and renamed to have a common format of your choosing. A FASTA file for targeted locus sequences can also be downloaded and included in the alignment.  See the specification of user-defined variables `mtGenPATTERN` and `GENBANKFASTA` below. 

Note that seaview is only used to convert from `fasta` to `nexus` format, so if you don't have it installed, you can manually convert the `fasta` to `nexus`. Also, while the pagan precompiled tar.gz does have mafft, it is not complete and you should use the complete mafft if you are aligning very long sequences, otherwise you will get an error when setting `LONGALIGNMENT=TRUE`

You can specify which positions to target (for specific loci and genes) to make alignments, including disjunct positions. For example, if you want to specify positions 1-10, then:

```bash
POSITIONS=1-10
```

If you want to align only positions 1-4 and 6-10, then:

```bash
POSITIONS=1-4,6-10
```

Refer to the mtGenomes annotation for the positions of particular loci of interest.

The other issue is which aligner to use.  I've tried `clustalw`, `clustalo`, `mafft`, and `pagan2`.  I've found `pagan2` to be superior in that it almost never needs to be aligned by eye to clean up mistakes.  The default behavior is to run `pagan2` for alignment.  However, if you run into problems, potentially due to sequences being very long, then you can try `mafft` with options set for very long sequences as follows:

```bash
LONGALIGNMENT=TRUE
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

bash radBARCODER.bash align $REF $bamPATTERN $THREADS $PREFIX $LOCUS $POSITIONS "$mtGenPATTERN" $LONGALIGNMENT $GENBANKFASTA
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

It is important to check the alignment by eye and edit as necessary or adjust upstream settings. I recommend [`seaview`](http://doua.prabi.fr/software/seaview) for this, but any alignment viewer will work. In the example data set, which has a lot of missing data, I did not have adjust the alignment at all, but mileage may vary.

#### 5. Lastly you can use `maximizeBP` to selectively cull your alignments from steps 4 or 7, either retaining more loci or more individuals, then goto step 6.

*Dependencies*: `R` (`seqinr`, `stringr`) 

This function will call `maximizeBP.R` which is included in the `radBARCODER` repo.  Make sure it is in your working directory

Set the `PCT` varable between 1 and 99, where it is the amount of allowable missing data. I recommend trying 10,25, and 50 to start with.  Histograms and culled alignments are output.  As the percent missing data goes down, the number of sequences retained also goes down, and the number of bp that are shared across all sequences goes up.

Update the following variable assignments and run `radBARCODER`:

```bash
FASTA=paganAlign_ALL_masked_aligned_clean_PfalcMitoGenome.fasta
PCT=99
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=95
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=90
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=85
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=80
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=75
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=70
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=50
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=30
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=25
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=20
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=15
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=10
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=5
bash radBARCODER.bash maximizeBP $FASTA $PCT
PCT=1
bash radBARCODER.bash maximizeBP $FASTA $PCT

```

Example successful output from 1 run looks like this (while some errors show up here, bute these are ok because they are not related to the `maximizeBP` function or are normal output that I have not suppressed (yet):

```
#########################################################################
Sat 05 Sep 2020 01:19:12 AM CDT RUNNING radBARCODER MAXIMIZEBP...
#########################################################################

Sat 05 Sep 2020 01:19:12 AM CDT VARIABLES READ IN:

the function that will be run FUNKTION=......maximizeBP
the reference genome used to map the reads REF=...........paganAlign_ALL_masked_aligned_clean_PfalcMitoGenome.fasta
the ls pattern shared by all bam files bamPATTERN=....99
the number of cpu cores for the task THREADS=.......
the characters added to every file created PREFIX=........
the name of the locus or loci LOCUS=.........
the nucleotide positions in the reference genome to consider POSITIONS=.....
the ls pattern shared by all mtGenomes that will be aligned mtGenPATTERN=..
the aligner that will be used LONGALIGNMENT=.
the GenBank sequences that should also be aligned GENBANK=.......

ls: cannot access '*99': No such file or directory

Sat 05 Sep 2020 01:19:12 AM CDT 1 SAMPLES BEING PROCESSED:



Sat 05 Sep 2020 01:19:12 AM CDT Be sure to check the fasta alignment by eye as small errors do occur
[1] "paganAlign_ALL_masked_aligned_clean_PfalcMitoGenome"
[2] "99"
null device
          1

#########################################################################
Sat 05 Sep 2020 01:19:13 AM CDT radBARCODER MAXIMIZEBP COMPLETED
#########################################################################
```

You can also check for proper output files to confirm success:

```
$ ls -tr | tail -n3
paganAlign_ALL_masked_aligned_clean_PfalcMitoGenome_99.fasta
missingCalls99.pdf
paganAlign_ALL_masked_aligned_clean_PfalcMitoGenome_99.nex
```


#### 6. Make network with `PopArt` 

[`PopArt`](https://github.com/jessicawleigh/popart-current), or your favorite network program, can now be used to create a network from the file.  `PopArt` automatically removes positions and sequences with poor coverage, so it's very convenient to apply to the file at this point.  [Precompiled, but outdated versions of PopArt](http://popart.otago.ac.nz/index.shtml)


#### 7. If you didn't have much luck comparing individuals in steps 1-5, you can make consensus sequences from groups of individuals and align those using `consensus` and then goto step 6

*Dependencies*: `R` (`seqinr`, `stringr`) 

Not vetted for mass consumption yet

This function will make consensus sequences for each sample category following the dDocent naming convention (`PopulationID_IndividualID`), but you need to specify the the population ids as described below.  The genesis of radBARCODER was trying to figure out what an unexpected population partition was, so it is also assumed that a subset of individuals will be identified at "nonTarget". A text file with one id (`PopulationID_IndividualID`) per line can be used for this as shown below.  The remaining individuals belonging to the majority or targeted taxon should be listed similarly in a separate file.
`
`cvgForCall` will determine the minimum read depth required to make a consensus base call. Base calling is performed by `consensusSeq.R` if you want to modify.

```bash
PREFIX=paganAlign_
LOCUS="PfalcMitoGenome"
THREADS=32
nontargetIDs=$(cat Pproctozystron.txt)
nontargetNAME=Ppr
targetIDs=$(cat Pfalcifer.txt)
targetNAME=Pfa
POPS=$(echo -e AtMk"\t"At"\t"Pk"\t"Kr"\t"St)
cvgForCall=1
bash radBARCODER.bash consensus "$nontargetIDs" "$targetIDs" $POPS $PREFIX $LOCUS $THREADS $cvgForCall $nontargetNAME $targetNAME
```

Intepreting errors: some error feedback is expected.  First, individuals that yielded no useful sequence are removed by `radBARCODER` and if they are listed as individuals from either the targeted or nontargeted taxon, they will trigger an error message, but will not affect the result.  





