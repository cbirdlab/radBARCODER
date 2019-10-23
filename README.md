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

You should `module load ddocent` on your hpc to get the required software for `bam2fasta`

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

Dependencies: ddocent pagan seaview 

```bash
#POSITIONS=60-550,6680-7020	#start and end positions of mtDNA fragment to excise, readable by cut -f 
#LOCUS="12S-COI"	#name of locus
#POSITIONS=40-200
POSITIONS=5665-5970
LOCUS="COI"#LOCUS="tRNA-Phe-12S"
#POSITIONS=10000-10500
#LOCUS="tRNA-Arg-ND4L-ND4"


#CUTOFFS=".Hspil.NC_023222"							#dDocent cutoffs used for reference genome
#PREFIX=Test	#prefix on files created
#THREADS=8   # number of cores to use
#mtGenPATTERN="reference.H*fasta"   #pattern match for fasta files with mito genomes to include in alignment
#GENBANKFASTA=""	#name of fasta file with additional sequences from genbank to include in alignment

CUTOFFS=".Hspil.NC_023222"
POSITIONS=1-17000
LOCUS="mtGenome"
PREFIX=Test_
THREADS=8
mtGenPATTERN="reference.H*fasta"
GENBANKFASTA=""
bash radBARCODER.bash align $CUTOFFS $THREADS $PREFIX $POSITIONS $LOCUS "$mtGenPATTERN" $GENBANKFASTA
```

#### 5. Make network with `PopArt` 

[`PopArt`](http://popart.otago.ac.nz/index.shtml), or your favorite network program, can now be used to create a network from the file.  `PopArt` automatically removes positions and sequences with poor coverage, so it's very convenient to apply to the file at this point.


#### 6. If you didn't have much luck comparing individuals in steps 1-5, you can make consensus sequences from groups of individuals and align those using `consensus` and then goto step 5


#### 7. Lastly you can use `maximizeBP` to selectively cull your alignments from steps 4 or 6, either retaining more loci or more individuals



