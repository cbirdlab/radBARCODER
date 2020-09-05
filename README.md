# radBARCODER

scripts to extract, align, and type mtDNA data from restriction site associated DNA sequenced on an [Illumina Machine](https://en.wikipedia.org/wiki/Illumina,_Inc.) with mitochondrial reference genomes of non-model species

---

## Preparing your [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files & a reference mitchondrial genome 

Follow these steps to make mtGenomes from each individual in your RAD data set.  We use the [dDocentHPC](https://github.com/cbirdlab/dDocentHPC) pipeline for processing RAD data in unix-based computers.  It is assumed that your FASTQ files are minimally processed (demultiplexed with no quality trimming) gzipped and have the following naming convention : 

```
# files must end with [FR].fq.gz
# only 1 underscore should occur, and it should delimit the population idenity and the individual identity.  
# every individual must have a different identity
Population_UniqueIndividualID.F.fq.gz
Population_UniqueIndividualID.R.fq.gz
```

It is also assumed that you have a fully assembled mitochondrial genome saved as a [FASTA](https://en.wikipedia.org/wiki/FASTA) file. You should name the reference mtGenome used for mapping sequence reads as follows:

```
# no more and no less than 3 periods should be used in the name and the * should be replaced with descriptive characters.
reference.*.*.fasta
```

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

## Quick Start

#### 1. Trim `fastq` files for mapping: [dDocentHPC trimFQmap](https://github.com/cbirdlab/dDocentHPC)

We will assume th

`config*` file settings:
```bash
32              Number of Processors (Auto, 1, 2, 3, ..., n threads) cbirdq=40 normal=20
230G    Maximum Memory (1G,2G,..., 256G)  G=gigabytes
----------trimFQ: Settings for Trimming FASTQ Files---------------------------------------------------------------
146             trimmomatic MINLEN (integer, mkREF only)                                                Drop the read if it is below a specified l$
75              trimmomatic MINLEN (integer, mkBAM only)                                                Drop the read if it is below a specified l$
20              trimmomatic LEADING:<quality> (integer, mkBAM only)                             Specifies the minimum quality required to keep a b$
15              trimmomatic TRAILING:<quality> (integer, mkREF only)                    Specifies the minimum quality required to keep a base.
20              trimmomatic TRAILING:<quality> (integer, mkBAM only)                    Specifies the minimum quality required to keep a base.
2               trimmomatic ILLUMINACLIP:<seed mismatches> (integer)                    specifies the maximum mismatch count which will still allo$
30              trimmomatic ILLUMINACLIP:<palindrome clip thresh> (integer)             specifies how accurate the match between the two 'adapter $
10              trimmomatic ILLUMINACLIP:<simple clip thresh> (integer)                 specifies how accurate the match between any adapter etc. $
20              trimmomatic SLIDINGWINDOW:<windowSize> (integer)                                specifies the number of bases to average across
20              trimmomatic SLIDINGWINDOW:<windowQuality> (integer)                             specifies the average quality required.
0               trimmomatic HEADCROP:<length> (integer, only Read1 for ezRAD)   The number of bases to remove from the start of the read. 0 for dd$
no              FixStacks (yes,no)                                                           Demultiplexing with stacks$
------------------------------------------------------------------------------------------------------------------
```

Run dDocentHPC to trim the reads as follows:

```bash
bash dDocentHPC.bash trimFQmap config.4.all
```

This will create a dir called `mkBAM` that is populated with trimmed `fq.gz` files.


#### 2. Map `fastq` to mtDNA genome using [dDocentHPC mkBAM](https://github.com/cbirdlab/dDocentHPC)

Move to the mkBAM dir

* obtain reference genome from [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
  * reference genome should be a `fasta` formatted file and can be composed of 1, several, or all loci in the mtGenome
  * name reference genome as follows: `reference.GenusSpecies.GenBankAccession.fasta` 
* set cutoff in the `dDocentHPC` `config*` file to *_GenusSpecies_*
* set cutoff2 i the `dDocentHPC` `config*` file to *_GenBankAccession_*

Here is an example of the `dDocentHPC` config file:

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

```bash
bash dDocentHPC.bash mkBAM config.4.all
bash dDocentHPC.bash fltrBAM config.4.all
```

This will create mildly filtered `RG.bam` files for each individual. These alignment maps are used in downstream processing.


#### 3. Create consensus sequences for each individual's reads mapped to the reference genome and mask areas with no coverage using `bam2fasta`

*Dependencies*: [`parallel`](https://www.gnu.org/software/parallel/) [`bedtools`](https://github.com/arq5x/bedtools2/releases) [`samtools`](https://www.htslib.org/)

From here forward, you'll be running the `radBARCODER.bash` script in `bash`.  It is highly likely that they can be completed on your laptop, but you will need the required software. 

You need to clone this repository to your computer:

```bash
git clone https://github.com/cbirdlab/radBARCODER.git   #html
```

Then you will need to copy all scripts in the repo ending in `bash` or `R` to your project directory:

```bash
#move into the repo you just cloned
cd PathTOradBARCODERrepo
cp *bash PathToProjectDir
cp *R PathToProjectDir
```

*Dependencies*: [`bedtools`](https://github.com/arq5x/bedtools2/releases) [`bcftools`](https://samtools.github.io/bcftools/bcftools.html) (fyi, both are required by ddocent, so you could `module load ddocent` on your hpc)

Update the following variable assignments and run `radBARCODER`:

```bash
REF=reference.Hspil.NC_023222.fasta  #Name of reference genome
bamPATTERN=.Hspil.NC_023222-RG.bam    #Pattern to id the bam files for each individual, this must be formatted as .CUTOFF1.CUTOFF2-RG.bam, where the CUTOFFs come from the dDocent config settings in step 2 above
THREADS=8                            #number of processors to use for parallel operations
bash radBARCODER.bash bam2fasta $REF $bamPATTERN $THREADS
```

This should result in a `vcf.gz` and a `masked_consensus.fasta` for every individual. Note that heterozygous positions are set to default to the reference allele. This behavior can be modified in `radBarcodder_functions.bash` at the line beginning with `bcftools consensus`. 

Hard coded strigencies are:

* depth of coverage >= 1

* base call phred quality >= 20

* mapping quality >= 30


#### 4. Select a portion of the genomes and `align` it across individuals

*Dependencies*: [`pagan`](http://wasabiapp.org/software/pagan/) [`mafft`](https://mafft.cbrc.jp/alignment/software/) [`seaview`](http://doua.prabi.fr/software/seaview) 

Note that seaview is only used to convert from `fasta` to `nexus` format, so if you don't have it installed, you can manually convert the `fasta` to `nexus`. Also, while the pagan precompiled tar.gz does have mafft, it is not complete and you should use the complete mafft if you are aligning very long sequences, otherwise you will get an error when setting `LONGALIGNMENT=TRUE`

You can specify which positions to target to make alignments, including disjunct positions. For example, if you want to specify positions 1-10, then:

```bash
POSITIONS=1-10
```

If you want positions 1-4 and 6-10, then:

```bash
POSITIONS=1-4,6-10
```

The other issue is which aligner to use.  I've tried `clustalw`, `clustalo`, `mafft`, and `pagan2`.  I've found `pagan2` to be superior in that it almost never needs to be aligned by eye to clean up mistakes.  The default behavior is to run `pagan2` for alignment.  However, if you run into problems, potentially due to sequences being very long, then you can try `mafft` with options set for very long sequences as follows:

```bash
LONGALIGNMENT=TRUE
```

Update the following variable assignments and run `radBARCODER`:

```bash
#REF=reference.Hspil.NC_023222.fasta  #Name of reference genome
#CUTOFFS=".Hspil.NC_023222"  #dDocent cutoffs used for reference genome
#PREFIX=Test  #prefix on files created
#THREADS=8  # number of cores to use
#mtGenPATTERN="*mtGenome.fasta"  #pattern match for fasta files with mito genomes to include in alignment
#GENBANKFASTA=""  #name of fasta file with additional sequences from genbank to include in alignment

REF=reference.Hspil.NC_023222.fasta  #Name of reference genome
bamPATTERN=.Hspil.NC_023222.bam
POSITIONS=40-200,5665-5970,10000-10500
LOCUS="tRNA-Phe-12S-COI-tRNA-Arg-ND4L-ND4"
PREFIX=Test_
THREADS=8
mtGenPATTERN="*mtGenome.fasta"
GENBANKFASTA=""

bash radBARCODER.bash align $REF $bamPATTERN $THREADS $PREFIX $LOCUS $POSITIONS "$mtGenPATTERN" $LONGALIGNMENT $GENBANKFASTA
```

It is important to check the alignment and edit as necessary. I recommend [`seaview`](http://doua.prabi.fr/software/seaview) for this, but any alignment viewer will work. 

#### 5. Make network with `PopArt` 

[`PopArt`](https://github.com/jessicawleigh/popart-current), or your favorite network program, can now be used to create a network from the file.  `PopArt` automatically removes positions and sequences with poor coverage, so it's very convenient to apply to the file at this point.  [Precompiled, but outdated versions of PopArt](http://popart.otago.ac.nz/index.shtml)


#### 6. If you didn't have much luck comparing individuals in steps 1-5, you can make consensus sequences from groups of individuals and align those using `consensus` and then goto step 5

*Dependencies*: `R` (`seqinr`, `stringr`) 

Not vetted for mass consumption yet

This function will make consensus sequences for each sample category following the dDocent naming convention (`PopulationID_IndividualID`), but you need to specify the the population ids as described below.  The genesis of radBARCODER was trying to figure out what an unexpected population partition was, so it is also assumed that a subset of individuals will be identified at "nonTarget". A text file with one id (`PopulationID_IndividualID`) per line can be used for this as shown below.  The remaining individuals belonging to the majority or targeted taxon should be listed similarly in a separate file.
`
`cvgForCall` will determine the minimum read depth required to make a consensus base call. Base calling is performed by `consensusSeq.R` if you want to modify.

```bash
PREFIX=concensusAlignment_
LOCUS="tRNA-Phe-12S-COI-tRNA-Arg-ND4L-ND4"
THREADS=32
nontargetIDs=$(cat RAD_OUTLIER_Pfalcifer_fish.txt)
targetIDs=$(cat RAD_NORMAL_Pfalcifer_fish.txt)
POPS=$(echo -e At"\t"Pk"\t"Kr"\t"St)
cvgForCall=1
bash radBARCODER.bash consensus $nontargetIDs $targetIDs $THREADS $PREFIX $LOCUS $POPS $cvgForCall
```

Intepreting errors: some error feedback is expected.  First, individuals that yielded no useful sequence are removed by `radBARCODER` and if they are listed as individuals from either the targeted or nontargeted taxon, they will trigger an error message, but will not affect the result.  


#### 7. Lastly you can use `maximizeBP` to selectively cull your alignments from steps 4 or 6, either retaining more loci or more individuals, then goto step 5.

*Dependencies*: `R` (`seqinr`, `stringr`) 

This function will call `maximizeBP.R` which is included in the `radBARCODER` repo.  Make sure it is in your working directory

Set the `PCT` varable between 1 and 99, where it is the amount of allowable missing data. I recommend trying 10,25, and 50 to start with.  Histograms and culled alignments are output.  As the percent missing data goes down, the number of sequences retained also goes down, and the number of bp that are shared across all sequences goes up.

Update the following variable assignments and run `radBARCODER`:

```bash
FASTA="Test_ALL_masked_aligned_clean_tRNA-Phe-12S-COI-tRNA-Arg-ND4L-ND4.fasta"
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

Despite the messages, if the files are created, then this worked as intended.


