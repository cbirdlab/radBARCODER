#!/bin/bash

#this script is meant to run the functions in findSpeciesID_functions.bash which enable the extraction of 
#mitochondrial data from ngs data, creation of consensus sequences, etc for the purposes of making
#haplotype networks and phylogenetic reconstructions

#source the functions
	source radBarcoder_functions.bash

#save commandline arguments to variables
	FUNKTION=$1
	REF=$2
	bamPATTERN=$3
	THREADS=$4
	PREFIX=$5
	LOCUS=$6
	POSITIONS=$7
	mtGenPATTERN=$8
	GENBANK=$9
	LONGALIGNMENT=$10

#user-defined variables
# if [ -z "$THREADS" ]; then THREADS=8 ; fi
# if [ -z "$CUTOFFS" ]; then CUTOFFS=".Hspil.NC_023222" ; fi	#dDocent cutoffs used for reference genome
# if [ -z "$PREFIX" ]; then PREFIX="" ; fi	#prefix on files created
# if [ -z "$PREFIX" ]; then PREFIX="" ; fi	#prefix on files created
# if [ -z "$PREFIX" ]; then PREFIX="" ; fi	#prefix on files created
# if [ -z "$PREFIX" ]; then PREFIX="" ; fi	#prefix on files created

#set other variables
#REF=reference${CUTOFFS}.fasta
#bamPATTERN=$CUTOFFS-RG
IDs=($(ls *$bamPATTERN | sed "s/$bamPATTERN//g" | grep -v '^cat'))
echo Samples being processed: 
echo ${IDs[@]}
echo ""
bamPATTERN=${bamPATTERN%.*}
LONGALIGNMENT=$(echo $LONGALIGNMENT | tr [a-z] [A-Z])

if [ "$FUNKTION" == "bam2fasta" ]; then
	#create fasta sequences of mtGenome for each individual using individually masked reference genomes
	echo ${IDs[@]} | tr " " "\n" | parallel -j $THREADS -k --no-notice "bam2fasta {} $bamPATTERN $REF"
	#concatenate the fasta mtGenome sequences from each fish
	cat *${bamPATTERN}_masked_consensus.fasta > all${bamPATTERN}_masked_consensus.fasta
elif [ "$FUNKTION" == "align" ]; then
	#get locus of choice from masked consensus seqs, mito genomes, and NCBI nucleotide seqs, clean and align
	alignLocusBySample $PREFIX $THREADS $bamPATTERN $POSITIONS $LOCUS "$mtGenPATTERN" $GENBANK $LONGALIGNMENT
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
