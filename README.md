# radBARCODER
scripts to extract, align, and type mtDNA data from restriction site associated DNA

---

## How it works

Follow these steps to make mtGenomes from each individual in your RAD data set.

#### 1. Trim `fastq` files for mapping: [dDocentHPC trimFQmap](https://github.com/cbirdlab/dDocentHPC)

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

This will create mildly filtered `RG.bam` files for each individual. These are your alignment maps are automatically used in downstream processing.


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

*Dependencies*: [`bedtools`](https://github.com/arq5x/bedtools2/releases) [`bcftools`](https://samtools.github.io/bcftools/bcftools.html) (fyi, both are required by ddocent, so you could module load ddocent on your hpc)

Update the following variable assignments and run `radBARCODER`:

```bash
REF=reference.Hspil.NC_023222.fasta  #Name of reference genome
bamPATTERN=.Hspil.NC_023222-RG.bam    #Pattern to id the bam files for each individual
THREADS=8                            #number of processors to use for parallel operations
bash radBARCODER.bash bam2fasta $REF $bamPATTERN $THREADS
```

This should result in a `vcf.gz` and a `masked_consensus.fasta` for every individual.  


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
#mtGenPATTERN="reference.H*fasta"  #pattern match for fasta files with mito genomes to include in alignment
#GENBANKFASTA=""  #name of fasta file with additional sequences from genbank to include in alignment

REF=reference.Hspil.NC_023222.fasta  #Name of reference genome
bamPATTERN=".Hspil.NC_023222.bam"
POSITIONS=40-200,5665-5970,10000-10500
LOCUS="tRNA-Phe-12S-COI-tRNA-Arg-ND4L-ND4"
PREFIX=Test_
THREADS=8
mtGenPATTERN="reference.H*fasta"
GENBANKFASTA=""

bash radBARCODER.bash align $REF $bamPATTERN $THREADS $PREFIX $LOCUS $POSITIONS "$mtGenPATTERN" $GENBANKFASTA
```

#### 5. Make network with `PopArt` 

[`PopArt`](https://github.com/jessicawleigh/popart-current), or your favorite network program, can now be used to create a network from the file.  `PopArt` automatically removes positions and sequences with poor coverage, so it's very convenient to apply to the file at this point.  [Precompiled, but outdated versions of PopArt](http://popart.otago.ac.nz/index.shtml)


#### 6. If you didn't have much luck comparing individuals in steps 1-5, you can make consensus sequences from groups of individuals and align those using `consensus` and then goto step 5

Not vetted for mass consumption yet

```bash
outIDs=$(cat RAD_OUTLIER_Pfalcifer_fish.txt)
normIDs=$(cat RAD_NORMAL_Pfalcifer_fish.txt)
SITES=$(echo -e At"\t"Pk"\t"Kr"\t"St)
cvgForCall=1
bash radBARCODER.bash consensus $outIDs $normIDs $THREADS $PREFIX $LOCUS $SITES $cvgForCall
```

#### 7. Lastly you can use `maximizeBP` to selectively cull your alignments from steps 4 or 6, either retaining more loci or more individuals, then goto step 5.

*Dependencies*: `R` (`seqinr`, `stringr`) 

This function will call `maximizeBP.R` which is included in the `radBARCODER` repo.  Make sure it is in your working directory

Set the `PCT` varable between 1 and 99, where it is the amount of allowable missing data. I recommend trying 10,25, and 50 to start with.  Histograms and culled alignments are output.  As the percent missing data goes down, the number of sequences retained also goes down, and the number of bp that are shared across all sequences goes up.

Update the following variable assignments and run `radBARCODER`:

```bash
FASTA="Test_ALL_masked_aligned_clean_tRNA-Phe-12S-COI-tRNA-Arg-ND4L-ND4.fasta"
PCT=10
bash radBARCODER.bash maximizeBP $FASTA $PCT
```



