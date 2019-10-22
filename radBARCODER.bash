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

#user-defined variables
if [ -z "$THREADS" ]; then THREADS=8 ; fi
if [ -z "$CUTOFFS" ]; then CUTOFFS=".Hspil.NC_023222" ; fi	#dDocent cutoffs used for reference genome


#set other variables
REF=reference${CUTOFFS}.fasta
bamPATTERN=$CUTOFFS-RG
IDs=($(ls *$bamPATTERN.bam | sed "s/$bamPATTERN\.bam//g"))


if FUNKTION=bam2fasta; then
	#create fasta sequences of mtGenome for each individual using individually masked reference genomes
	echo ${IDs[@]} | tr " " "\n" | parallel -j $THREADS -k --no-notice "bam2fasta {} $bamPATTERN $REF"
	#concatenate the fasta mtGenome sequences from each fish
	cat *${bamPATTERN}_masked_consensus.fasta > all${bamPATTERN}_masked_consensus.fasta
fi
