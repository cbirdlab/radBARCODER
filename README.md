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

Run dDocentHPC as follows:
```bash
bash dDocentHPC.bash trimFQmap config.4.all
```

#### 2. Map `fastq` to mtDNA genome using [dDocentHPC mkBAM](https://github.com/cbirdlab/dDocentHPC)
  * obtain reference genome from [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
    * reference genome should be a `fasta` formatted file and can be composed of 1, several, or all loci in the mtGenome
    * name reference genome as follows: `reference.GenusSpecies.GenBankAccession.fasta` 
  * set cutoff in the `dDocentHPC` `config*` file to *_GenusSpecies_*
  * set cutoff2 i the `dDocentHPC` `config*` file to *_GenBankAccession_*

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

#### 3. Create consensus sequences for each individual's reads mapped to the reference genome and mask areas with no coverage using `bam2fasta`

dependencies: bedtools bcftools (fyi, both are required by ddocent)

```bash
CUTOFFS=".Hspil.NC_023222"							#dDocent cutoffs used for reference genome
THREADS=8
#bamPATTERN=$CUTOFFS-RG           #search pattern for bam files
#REF=reference${CUTOFFS}.fasta
#IDs=($(ls *$bamPATTERN.bam | sed "s/$bamPATTERN\.bam//g"))

bash radBARCODER.bash bam2fasta $CUTOFFS $THREADS

```

This should result in a `vcf.gz` and a `masked_consensus.fasta` for every individual.  

#### 4. Select a portion of the genomes and `align` it across individuals

Dependencies: pagan mafft seaview 

You can specify which positions to target to make alignments, including disjunct positions. For example, if you want to specify positions 1-10, then:

```bash
POSITIONS=1-10
```

If you want positions 1-4 and 6-10, then:

```bash
POSITIONS=1-4,6-10
```

The other issue is which aligner to use.  I've tried clustalw, clustalo, mafft, and pagan2.  I've found pagan2 to be superior in that it almost never needs to be aligned by eye to clean up mistakes.  The default behavior is to run pagan2 for alignment.  However, if you run into problems, potentially due to sequences being very long, then you can try mafft with options set for very long sequences as follows:

```bash
LONGALIGNMENT=TRUE
```

The following are the variables to set and command to run the alignment.

```bash
POSITIONS=60-550,6680-7020	#start and end positions of mtDNA fragment to excise, readable by cut -f 
LOCUS="12S-COI"	#name of locus
POSITIONS=40-200
LOCUS="tRNA-Phe-12S"
POSITIONS=5665-5970
LOCUS="COI"
POSITIONS=10000-10500
LOCUS="tRNA-Arg-ND4L-ND4"
POSITIONS=1-17000
LOCUS="mtGenome"


#CUTOFFS=".Hspil.NC_023222"							#dDocent cutoffs used for reference genome
#PREFIX=Test	#prefix on files created
#THREADS=8   # number of cores to use
#mtGenPATTERN="reference.H*fasta"   #pattern match for fasta files with mito genomes to include in alignment
#GENBANKFASTA=""	#name of fasta file with additional sequences from genbank to include in alignment

CUTOFFS=".Hspil.NC_023222"
POSITIONS=40-200,5665-5970,10000-10500
LOCUS="tRNA-Phe-12S-COI-tRNA-Arg-ND4L-ND4"
PREFIX=Test_
THREADS=8
mtGenPATTERN="reference.H*fasta"
GENBANKFASTA=""
LONGALIGNMENT=FALSE
bash radBARCODER.bash align $CUTOFFS $THREADS $PREFIX $LOCUS $POSITIONS "$mtGenPATTERN" $GENBANKFASTA
```

#### 5. Make network with `PopArt` 

[`PopArt`](https://github.com/jessicawleigh/popart-current), or your favorite network program, can now be used to create a network from the file.  `PopArt` automatically removes positions and sequences with poor coverage, so it's very convenient to apply to the file at this point.  [Precompiled, but outdated versions of PopArt](http://popart.otago.ac.nz/index.shtml)


#### 6. If you didn't have much luck comparing individuals in steps 1-5, you can make consensus sequences from groups of individuals and align those using `consensus` and then goto step 5


#### 7. Lastly you can use `maximizeBP` to selectively cull your alignments from steps 4 or 6, either retaining more loci or more individuals

dependencies: R (seqinr, stringr)

Set the PCT varable between 1 and 99, where it is the amount of allowable missing data. I recommend trying 10,25, and 50 to start with.  Histograms and culled alignments are output.  As the percent missing data goes down, the number of sequences retained also goes down, and the number of bp that are shared across all sequences goes up.

```bash
FASTA="Test_ALL_masked_aligned_clean_tRNA-Phe-12S-COI-tRNA-Arg-ND4L-ND4.fasta"
PCT=10
bash radBARCODER.bash maximizeBP $FASTA $PCT
```



