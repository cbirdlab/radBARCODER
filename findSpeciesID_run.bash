#!/bin/bash

#this script is meant to run the functions in findSpeciesID_functions.bash which enable the extraction of 
#mitochondrial data from ngs data, creation of consensus sequences, etc for the purposes of making
#haplotype networks and phylogenetic reconstructions

#source the functions
	source findSpeciesID_functions.bash

#user-defined variables
	THREADS=7
	CUTOFFS=".48.48"							#dDocent cutoffs used for reference genome
	FILTER=4									#filter reads with this samtools flag
	PREFIX=fish									#prefix on files created
	
	LOCUS="mtDNA2_1"									#name of locus
	POSITIONS=1-16000							#start and end positions of mtDNA fragment to excise, readable by cut -f 
	GENBANK=""
	
	# LOCUS="mtDNA"									#name of locus
	# POSITIONS=1-7000							#start and end positions of mtDNA fragment to excise, readable by cut -f 
	# GENBANK=""

	# #high rad coverage
	# LOCUS="tRNAPhe-12S"									#name of locus
	# POSITIONS=60-550						#start and end positions of mtDNA fragment to excise, readable by cut -f 
	# GENBANK=""

	# # high rad coverage
	# LOCUS="COI"									#name of locus
	# POSITIONS=6680-7020						#start and end positions of mtDNA fragment to excise, readable by cut -f 
	# GENBANK=""

	mtGENOMEs="*Puntioplites*mtGenome.fasta"	#ls search string to return mitoGenomes, names should end with _mtGenome.fasta
	OUTLIERS=RAD_OUTLIER_Pfalcifer_fish.txt		#name of file containing fish in outlier group
	NORMALS=RAD_NORMAL_Pfalcifer_fish.txt		#name of file that will be made containing fish in normal group
	SITES=("At" "Kr" "Pk" "St")					#used to make consensus for each site, seq names begin with these letters
	cvgForCall=3								#number of individuals with ACTG at pos to induce consensus call, majority
	
#set other variables
	REF=reference${CUTOFFS}.fasta
	midFILE=$CUTOFFS-FLTR_F$FILTER
	IDs=($(ls *_Pfa*$midFILE.bam | sed "s/$midFILE\.bam//g"))
	outIDs=($(cat $OUTLIERS))
	normIDs=($(echo ${IDs[@]} | tr " " "\n" | grep -v -f $OUTLIERS))
	mtGenNames=($(ls $mtGENOMEs))
	
#make files used by functions
	echo ${normIDs[@]} | tr " " "\n" > $NORMALS
	

#create fasta sequences of mtGenome for each individual using individually masked reference genomes
	echo ${IDs[@]} | tr " " "\n" | parallel -j $THREADS -k --no-notice "bam2fasta {} $midFILE $REF"
	#concatenate the fasta mtGenome sequences from each fish
		cat *${midFILE}_masked_consensus.fasta > all${midFILE}_masked_consensus.fasta

#get locus of choice from masked consensus seqs, mito genomes, and NCBI nucleotide seqs, clean and align
	alignLocusBySample $PREFIX $THREADS $midFILE $POSITIONS $LOCUS $GENBANK IDs

#make consensus sequences from aligned fasta files
	mkConsensusFasta outIDs normIDs SITES $PREFIX $LOCUS $THREADS $cvgForCall



#merge all bam files
#catBAM $CUTOFFS 

#use whole cat file to get consensus seqs for each sample then extract a locus and align
#getConsensusBySample $CUTOFFS $FILTER $PREFIX 
#alignBySample $CUTOFFS $FILTER $PREFIX $POSITIONS $LOCUS $mtGENOMEs

#merge all bam files by outlier group or normal group
#grpCatBAM $CUTOFFS $OUTLIERS

#use whole cat file to get consensus seqs for each group then extract a locus and align
#getConsensusByGroup $CUTOFFS $FILTER $PREFIX 
#alignByGroup $CUTOFFS $FILTER $PREFIX $POSITIONS $LOCUS $mtGENOMEs
