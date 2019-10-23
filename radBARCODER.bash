#!/bin/bash

#this script is meant to run the functions in findSpeciesID_functions.bash which enable the extraction of 
#mitochondrial data from ngs data, creation of consensus sequences, etc for the purposes of making
#haplotype networks and phylogenetic reconstructions

#source the functions
	source radBarcoder_functions.bash

#save commandline arguments to variables
	FUNKTION=$1
	CUTOFFS=$2
	THREADS=$3
	PREFIX=$4
	LOCUS=$5
	POSITIONS=$6
	mtGenPATTERN=$7
	GENBANK=$8
	LONGALIGNMENT=$9

#user-defined variables
# if [ -z "$THREADS" ]; then THREADS=8 ; fi
# if [ -z "$CUTOFFS" ]; then CUTOFFS=".Hspil.NC_023222" ; fi	#dDocent cutoffs used for reference genome
# if [ -z "$PREFIX" ]; then PREFIX="" ; fi	#prefix on files created
# if [ -z "$PREFIX" ]; then PREFIX="" ; fi	#prefix on files created
# if [ -z "$PREFIX" ]; then PREFIX="" ; fi	#prefix on files created
# if [ -z "$PREFIX" ]; then PREFIX="" ; fi	#prefix on files created

#set other variables
REF=reference${CUTOFFS}.fasta
bamPATTERN=$CUTOFFS-RG
IDs=($(ls *$bamPATTERN.bam | sed "s/$bamPATTERN\.bam//g"))
LONGALIGNMENT=$(echo $LONGALIGNMENT | tr [a-z] [A-Z])

if [ "$FUNKTION" == "bam2fasta" ]; then
	#create fasta sequences of mtGenome for each individual using individually masked reference genomes
	echo ${IDs} | tr " " "\n" | parallel -j $THREADS -k --no-notice "bam2fasta {} $bamPATTERN $REF"
	#concatenate the fasta mtGenome sequences from each fish
	cat *${bamPATTERN}_masked_consensus.fasta > all${bamPATTERN}_masked_consensus.fasta
elif [ "$FUNKTION" == "align" ]; then
	#get locus of choice from masked consensus seqs, mito genomes, and NCBI nucleotide seqs, clean and align
	alignLocusBySample $PREFIX $THREADS $bamPATTERN $POSITIONS $LOCUS "$mtGenPATTERN" $GENBANK
elif [ "$FUNKTION" == "consensus" ]; then
	#make consensus sequences from aligned fasta files
		mkConsensusFasta outIDs normIDs SITES $PREFIX $LOCUS $THREADS $cvgForCall
elif [ "$FUNKTION" == "maximizeBP" ]; then
	#maximize the number of bp retained at the expense of retaining individuals
		FASTA=$2
		PCT=$3
		echo; echo `date` "Be sure to check the fasta alignment by eye as small errors do occur"
		maximizeBP $FASTA $PCT
fi