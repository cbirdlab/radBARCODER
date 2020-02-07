#!/bin/bash

#this script contains functions to convert bam files to fasta alignments for the purposes of DNA barcoding

#default user defined variables
# THREADS=8
# CUTOFFS=".Hspil.NC_023222"							#dDocent cutoffs used for reference genome
#FILTER=4									#filter reads with this samtools flag
# PREFIX=Hspilopterus									#prefix on files created
# POSITIONS=5683-5951							#start and end positions of mtDNA fragment to excise
# LOCUS=begCOI								#name of locus
# mtGENOMEs="reference.fasta"	#ls search string to return mitoGenomes, names should end with _mtGenome.fasta
# OUTLIERS=RAD_OUTLIER_Pfalcifer_fish.txt		#name of file containing fish in outlier group
pctMissCall=50

#automatic variables
# REF=reference${CUTOFFS}.fasta
BAMLIST=bamlist${CUTOFFS}.list
catFILE1=cat${CUTOFFS}-FLTR_F4.bam 
bamLIST2=bamlist${CUTOFFS}.list2
catFILE2=cat${CUTOFFS}-FLTR_F$FILTER.bam 


#
ntrlvdFasta2Fasta(){
	#convert an interleaved fasta to a normal fasta
	if read -t 0; then
	INFILE=$1
	OUTFILE=$2
	SORT=$3
	if [ "$SORT" == "TRUE" ]; then
		cat $INFILE | tr -d "\r" | sed 's/\(^>.*\)$/\1\t/' | tr -d "\n" | sed -e $'s/>/\\n>/g' | tail -n+2 | sort | tr "\t" "\n" > $OUTFILE
	else
		cat $INFILE | tr -d "\r" | sed 's/\(^>.*\)$/\1\t/' | tr -d "\n" | sed -e $'s/>/\\n>/g' | tail -n+2 | tr "\t" "\n" > $OUTFILE
	fi

	else

	INFILE=/dev/stdin
	SORT=$1
	if [ "$SORT" == "TRUE" ]; then
		cat $INFILE | tr -d "\r" | sed 's/\(^>.*\)$/\1\t/' | tr -d "\n" | sed -e $'s/>/\\n>/g' | tail -n+2 | sort | tr "\t" "\n" 
	else
		cat $INFILE | tr -d "\r" | sed 's/\(^>.*\)$/\1\t/' | tr -d "\n" | sed -e $'s/>/\\n>/g' | tail -n+2 | tr "\t" "\n" 
	fi

	fi
}
export -f ntrlvdFasta2Fasta


#function to create consensus sequences using individually masked reference genomes
	#superior to using a single unmasked reference genome
	#designed to be run in parallel
bam2fasta(){
	#assign arguments to variables
		local ID=$1							#*_Pfa*
		local ID2=$1$2						#*_Pfa*$CUTOFFS-FLTR_F4
		local REF=$3						#unmasked reference genome
	#calculate coverage at each base from bam file and return a bed file with positions that have no coverage
		bedtools genomecov -ibam $ID2.bam -bga | grep -P '\t0$' > $ID2.bed 
	#create masked fasta for each individual
		bedtools maskfasta -fi $REF -fo ${ID2}_masked_ref.fasta -bed $ID2.bed 
	#make vcf files
		bcftools mpileup --threads 7 -d 30000 -q 30 -Q 20 -A -O z -o ${ID2}_masked_pile.vcf.gz -f ${ID2}_masked_ref.fasta ${ID2}.bam
	#call genotypes in vcf, force ploidy=haploid
		bcftools call --threads 7 -m --ploidy 1 -O z -o ${ID2}_masked_calls.vcf.gz ${ID2}_masked_pile.vcf.gz
	#combine snps and indels into 1 multiallelic call
		bcftools norm -f ${ID2}_masked_ref.fasta -m +any -O z -o ${ID2}_masked_calls_normalized.vcf.gz ${ID2}_masked_calls.vcf.gz
		tabix ${ID2}_masked_calls_normalized.vcf.gz
	#generate 1 consensus for each fishSample
		bcftools consensus -s $ID -M 'N' -f ${ID2}_masked_ref.fasta ${ID2}_masked_calls_normalized.vcf.gz -o ${ID2}_masked_consensus.fasta
	#add ID to name of consensus
		sed -i "s/^>/>${ID}_/g" ${ID2}_masked_consensus.fasta
}
export -f bam2fasta

#function to get locus from masked consensus sequences, mito genomes, and NCBI nucleotide records, clean and align
alignLocusBySample(){
	#assign arguments to variables
		local PREFIX=$1
		local THREADS=$2
		local midFILE=$3
		local POSITIONS=$4
		local LOCUS=$5
		local mtGenPATTERN="$6"		#cant figure this one out
		local LONGALIGNMENT=$7
		local GENBANK=$8

		# local name=$7[@]
		# local IDs=("${!name}")

	#get locus from all individuals and align, $POSITIONS determines the locus
		echo ${IDs[@]} | tr " " "\n" | parallel -j $THREADS -k --no-notice "echo \>{} && tail -n +2 {}${midFILE}_masked_consensus.fasta | tr '\n' '\t' | sed 's/\t//g' | sed 's/ */\t/g' | cut -f $POSITIONS | sed 's/\t//g' " > ${PREFIX}RAD_masked_$LOCUS.fasta
	#remove individuals with no sequences or all N
		sed 's/^NN*$//g' ${PREFIX}RAD_masked_$LOCUS.fasta | grep -P -B 1 '^[ACTGN@]+$' | grep -v '\-\-' > ${PREFIX}RAD_masked_${LOCUS}_clean.fasta
	#make fasta from mtGenome sequences
		if [ ! -z "$mtGenPATTERN" ]; then
			ls $mtGenPATTERN | sed 's/.fasta//g' | parallel -j $THREADS -k --no-notice "echo \>{} && tail -n +2 {}.fasta | tr '\n' '\t' | sed 's/\t//g' | sed 's/ */\t/g' | cut -f $POSITIONS | sed 's/\t//g' " > ${PREFIX}MtGenomes_$LOCUS.fasta
		fi
	#Make fasta from other NCBI nucleotide sequences
		#need to code, right now, just download as 1 big fasta
	#Combine fastas, make sure to add genbank nucleotide recs as neccessary
		#echo $GENBANK
		cat ${PREFIX}MtGenomes_$LOCUS.fasta ${PREFIX}RAD_masked_${LOCUS}_clean.fasta $GENBANK  > ${PREFIX}ALL_masked_$LOCUS.fasta 
	#Align all_*.fasta
		#clustalw -infile=${PREFIX}ALL_masked_$LOCUS.fasta -align -type=DNA -output=NEXUS -outfile=${PREFIX}ALL_masked_aligned_$LOCUS.nex
		#clustalo -infile=${PREFIX}ALL_masked_$LOCUS.fasta -t DNA --outfmt=fa -outfile=${PREFIX}ALL_masked_aligned_$LOCUS.fasta --threads $THREADS 
		#mafft --thread $THREADS ${PREFIX}ALL_masked_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_$LOCUS.fasta
		#mafft --thread $THREADS --ep 0.123 ${PREFIX}ALL_masked_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_$LOCUS.fasta
		echo LONGALIGNMENT=$LONGALIGNMENT
		if [ "$LONGALIGNMENT" == "TRUE" ] || [ "$LONGALIGNMENT" == "T" ]; then
			mafft --thread $THREADS --globalpair --maxiterate 1000 ${PREFIX}ALL_masked_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_$LOCUS.fasta
		else 
			pagan2 -s ${PREFIX}ALL_masked_$LOCUS.fasta -o ${PREFIX}ALL_masked_aligned_$LOCUS --threads $THREADS
			mv ${PREFIX}ALL_masked_aligned_$LOCUS.fas ${PREFIX}ALL_masked_aligned_${LOCUS}.fasta
		fi
	#clean out sites with all gaps and 'n's, make sure to skip the lines with names
		sed -e '/>/!s/[nN]/\-/g' ${PREFIX}ALL_masked_aligned_$LOCUS.fasta | \
		seaview -convert -output_format fasta -o ${PREFIX}ALL_masked_aligned_clean_$LOCUS.fasta -del_gap_only_sites -
		seaview -convert -output_format nexus -o ${PREFIX}ALL_masked_aligned_clean_$LOCUS.nex ${PREFIX}ALL_masked_aligned_clean_$LOCUS.fasta
}
export -f alignLocusBySample

#function to make consensus sequences from aligned fasta files
mkConsensusFasta(){
	#assign arguments to variables
		name=$1[@]
		local outIDs=("${!name}")
		name=$2[@]
		local normIDs=("${!name}")
		name=$3[@]
		local SITES=("${!name}")
		local PREFIX=$4
		local LOCUS=$5
		local THREADS=$6
		local cvgForCall=$7
	
	#split up fasta for making consensus seqs and name files according to seq id
		#split up fasta by sequence and make list of file names
			csplit --quiet --digits=4 --prefix=split${PREFIX}_masked_aligned_$LOCUS.fasta ${PREFIX}ALL_masked_aligned_$LOCUS.fasta "/^>/+0" "{*}"
			rm split${PREFIX}_masked_aligned_$LOCUS.fasta0000
			local fileNAMES=($(ls split${PREFIX}_masked_aligned_$LOCUS.fasta*))
			local seqNAMES=($(parallel -j $THREADS -k --no-notice "grep '^>' {} | sed 's/>//g' " ::: ${fileNAMES[@]}))
		#rename the files by sequence name
			parallel -j $THREADS -k --no-notice --link "mv {1} ${PREFIX}{2}_masked_aligned_$LOCUS.fasta" ::: ${fileNAMES[@]} ::: ${seqNAMES[@]} 
			parallel -j $THREADS -k --no-notice --link "echo {1} {2}" ::: ${fileNAMES[@]} ::: ${seqNAMES[@]} 

	#make consensus sequence from outlier fish
		parallel -j $THREADS -k --no-notice "cat ${PREFIX}{}_masked_aligned_$LOCUS.fasta" ::: ${outIDs[@]} | grep -v '^cat:' | grep -v '^cat:' > ${PREFIX}OUTLIERS_masked_aligned_$LOCUS.fasta
		Rscript consensusSeq.R ${PREFIX}OUTLIERS_masked_aligned_$LOCUS.fasta ${PREFIX}OUTLIERS_masked_aligned_consensus_$LOCUS.fasta $cvgForCall Outlier_Consensus

	#make consensus sequence from normal fish
		parallel -j $THREADS -k --no-notice "cat ${PREFIX}{}_masked_aligned_$LOCUS.fasta" ::: ${normIDs[@]} | grep -v '^cat:' | grep -v '^cat:' > ${PREFIX}NORMALS_masked_aligned_$LOCUS.fasta
		Rscript consensusSeq.R ${PREFIX}NORMALS_masked_aligned_$LOCUS.fasta ${PREFIX}NORMALS_masked_aligned_consensus_$LOCUS.fasta $cvgForCall Normal_Consensus

	#make consensus sequence from normal fish for each sample location
		ncbiNAMES=($(echo ${seqNAMES[@]} | tr " " "\n" ))
		for i in ${SITES[@]}; do
			siteIDs=($(echo ${normIDs[@]} | tr " " "\n" | grep "^$i"))
			parallel -j $THREADS -k --no-notice "cat ${PREFIX}{}_masked_aligned_$LOCUS.fasta" ::: ${siteIDs[@]} | grep -v '^cat:' | grep -v '^cat:' > ${PREFIX}${i}_masked_aligned_$LOCUS.fasta
			Rscript consensusSeq.R ${PREFIX}${i}_masked_aligned_$LOCUS.fasta ${PREFIX}${i}_masked_aligned_consensus_$LOCUS.fasta $cvgForCall ${i}_Consensus

			#get list of NCBI seqs
				ncbiNAMES=($(echo ${ncbiNAMES[@]} | tr " " "\n" | grep -v "$i"))
		done
	
	#get ncbi seqs
		parallel -j $THREADS -k --no-notice "cat ${PREFIX}{}_masked_aligned_$LOCUS.fasta" ::: ${ncbiNAMES[@]} > ${PREFIX}NCBI_masked_aligned_$LOCUS.fasta
	
	#concatenate consensus sequences into 1 file
		cat ${PREFIX}NCBI_masked_aligned_$LOCUS.fasta ${PREFIX}*_masked_aligned_consensus_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_consensi_$LOCUS.fasta
	#mafft makes a bunch of sites with all indels due to a few poorly aligned sequences, remove gap only sites
		seaview -convert -output_format fasta -o ${PREFIX}ALL_masked_aligned_clean_consensi_$LOCUS.fasta -del_gap_only_sites ${PREFIX}ALL_masked_aligned_consensi_$LOCUS.fasta
	#convert fasta to nexus non interleaved
		seaview -convert -output_format nexus -o ${PREFIX}ALL_masked_aligned_clean_consensi_$LOCUS.nex ${PREFIX}ALL_masked_aligned_clean_consensi_$LOCUS.fasta
}

#function to maximize the number of bp retained at the expense of retaining individuals
maximizeBP(){
	#read arguments into variables
	local INFILE=$1
	local pctMissCall=$2
	Rscript maximizeBP.R ${INFILE%.*} $pctMissCall
	seaview -convert -output_format nexus -o ${INFILE%.*}_$pctMissCall.nex ${INFILE%.*}_$pctMissCall.fasta
}



