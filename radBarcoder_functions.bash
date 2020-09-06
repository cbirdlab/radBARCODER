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


# function to create consensus sequences using individually masked reference genomes
	# superior to using a single unmasked reference genome
	# designed to be run in parallel
bam2fasta(){
	# assign arguments to variables
		local ID=$1							#*_Pfa*
		local ID2=$1$2						#*_Pfa*$CUTOFFS-FLTR_F4
		local REF=$3						#unmasked reference genome
	# calculate coverage at each base from bam file and return a bed file with positions that have no coverage
		bedtools genomecov -ibam $ID2.bam -bga | grep -P '\t0$' > $ID2.bed
	# create masked fasta for each individual
		bedtools maskfasta -fi $REF -fo ${ID2}_masked_ref.fasta -bed $ID2.bed
	# make vcf files
		bcftools mpileup --threads 1 -d 30000 -q 30 -Q 20 -A -O z -o ${ID2}_masked_pile.vcf.gz -f ${ID2}_masked_ref.fasta ${ID2}.bam
	# call genotypes in vcf, force ploidy=haploid
		bcftools call --threads 1 -m --ploidy 1 -O z -o ${ID2}_masked_calls.vcf.gz ${ID2}_masked_pile.vcf.gz
	# combine snps and indels into 1 multiallelic call
		bcftools norm -f ${ID2}_masked_ref.fasta -m +any -O z -o ${ID2}_masked_calls_normalized.vcf.gz ${ID2}_masked_calls.vcf.gz
		tabix ${ID2}_masked_calls_normalized.vcf.gz
	# generate 1 consensus for each fishSample
		bcftools consensus -s $ID -M 'N' -f ${ID2}_masked_ref.fasta ${ID2}_masked_calls_normalized.vcf.gz -o ${ID2}_masked_consensus.fasta
	# add ID to name of consensus
		sed -i "s/^>/>${ID}_/g" ${ID2}_masked_consensus.fasta
	# clean up files
		mkdir bam2fasta_out
		mv $ID2.bed bam2fasta_out
		mv ${ID2}_masked_ref.fasta* bam2fasta_out
		mv ${ID2}_masked_calls_normalized.vcf.gz* bam2fasta_out
		mv ${ID2}_masked_calls.vcf.gz* bam2fasta_out
		mv ${ID2}_masked_pile.vcf.gz* bam2fasta_out
		
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
		local REF=$9

		# local name=$7[@]
		# local IDs=("${!name}")

	#get locus from all individuals and align, $POSITIONS determines the locus
	echo ""; echo `date` EXTRACTING POSITIONS $POSITIONS FOR ALIGNMENT...
		echo ${IDs[@]} | tr " " "\n" | parallel -j $THREADS -k --no-notice "echo \>{} && tail -n +2 {}${midFILE}_masked_consensus.fasta | tr '\n' '\t' | sed 's/\t//g' | sed 's/ */\t/g' | cut -f $POSITIONS | sed 's/\t//g' " > ${PREFIX}RAD_masked_$LOCUS.fasta

	#remove individuals with no sequences or all N
	echo ""; echo `date` REMOVING INDIVIDUALS WITH NO NUCLEOTIDES CALLED...
		sed 's/^NN*$//g' ${PREFIX}RAD_masked_$LOCUS.fasta | grep -P -B 1 '^[ACTGN@]+$' | grep -v '\-\-' > ${PREFIX}RAD_masked_${LOCUS}_clean.fasta

	#make fasta from mtGenome sequences
	if [ ! -z "$mtGenPATTERN" ]; then
		echo ""; echo `date` GATHERING ALL mtGENOMEs WITH PATTERN=$mtGenPATTERN AND EXTRACTING POSITIONS $POSITIONS
		ls $mtGenPATTERN | sed 's/.fasta//g' | parallel -j $THREADS -k --no-notice "echo \>{} && tail -n +2 {}.fasta | tr '\n' '\t' | sed 's/\t//g' | sed 's/ */\t/g' | cut -f $POSITIONS | sed 's/\t//g' " > ${PREFIX}MtGenomes_$LOCUS.fasta
		LinesInMtGenomes=$(wc -l ${PREFIX}MtGenomes_$LOCUS.fasta | cut -d" " -f1)
	else
		echo ""; echo `date` NO mtGENOMES SPECIFIED, CONTINUING WITHOUT THEM...
	fi
	#Make fasta from other NCBI nucleotide sequences
		#need to code, right now, just download as 1 big fasta

	#Combine fastas, make sure to add genbank nucleotide recs as neccessary
		if [ ! -z "$GENBANK" ]; then
			echo ""; echo `date` ADDING $GENBANK SEQUENCES FROM GENBANK...
		else
			echo ""; echo `date` NO ADDITIONAL GENBANK SEQUENCES SPECIFIED, CONTINUING WITHOUT THEM...
		fi
		echo ""; echo `date` CONCATENATING FASTAs...
		cat ${PREFIX}MtGenomes_$LOCUS.fasta ${PREFIX}RAD_masked_${LOCUS}_clean.fasta $GENBANK  > ${PREFIX}ALL_masked_$LOCUS.fasta
	#Align all_*.fasta
	echo ""; echo `date` ALIGNING SEQUENCES WITH $(if [ "$LONGALIGNMENT" == "TRUE" ]; then echo -n MAFFT; else echo -n pagan2; fi)...
		#clustalw -infile=${PREFIX}ALL_masked_$LOCUS.fasta -align -type=DNA -output=NEXUS -outfile=${PREFIX}ALL_masked_aligned_$LOCUS.nex
		#clustalo -infile=${PREFIX}ALL_masked_$LOCUS.fasta -t DNA --outfmt=fa -outfile=${PREFIX}ALL_masked_aligned_$LOCUS.fasta --threads $THREADS
		#mafft --thread $THREADS ${PREFIX}ALL_masked_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_$LOCUS.fasta
		#mafft --thread $THREADS --ep 0.123 ${PREFIX}ALL_masked_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_$LOCUS.fasta
		echo ""
		if [ "$LONGALIGNMENT" == "TRUE" ] || [ "$LONGALIGNMENT" == "T" ]; then
			mafft --thread $THREADS --globalpair --maxiterate 1000 ${PREFIX}ALL_masked_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_$LOCUS.fasta
		else
			# my first try at making this work
			# pagan2 -s ${PREFIX}ALL_masked_$LOCUS.fasta -o ${PREFIX}ALL_masked_aligned_$LOCUS --threads $THREADS

			# second method from author of pagan, getting run out of memory error
			#sed 's/N\{20,\}/N/g' ${PREFIX}RAD_masked_${LOCUS}_clean.fasta > ${PREFIX}RAD_masked_${LOCUS}_clean_2.fasta
			#pagan2 -s ${PREFIX}MtGenomes_$LOCUS.fasta -o ref
			#pagan2 -a ref.fas -r ref.tre -q ${PREFIX}RAD_masked_${LOCUS}_clean_2.fasta -o ${PREFIX}ALL_masked_aligned_$LOCUS --pileup --no-terminal-edges

			# third method from author of pagan
			#sed -i 's/N\{20,\}/N/g' ${PREFIX}ALL_masked_${LOCUS}.fasta
			#pagan2 -q ${PREFIX}ALL_masked_$LOCUS.fasta -o ${PREFIX}ALL_masked_aligned_$LOCUS --pileup --no-terminal-edges

			# my solution for pagan alignment: align 1 individual at a time to the mulitiple ref genomes
			echo `date` ALIGNING RAD DATA TO MITOGENOMES...
			sed 's/N\{20,\}/N/g' ${PREFIX}RAD_masked_${LOCUS}_clean.fasta > ${PREFIX}RAD_masked_${LOCUS}_clean_2.fasta
			IndivNames=($(cat ${PREFIX}RAD_masked_${LOCUS}_clean_2.fasta | paste - - | sed 's/^>//' | cut -f1))
			IndivSeqs=($(cat ${PREFIX}RAD_masked_${LOCUS}_clean_2.fasta | paste - - | sed 's/^>//' | cut -f2))
			parallel --no-notice --link -j $THREADS "printf '>{1}\n{2}\n' > PrePaganAligned_{1}.fasta " ::: ${IndivNames[@]} ::: ${IndivSeqs[@]}
			pagan2 -s ${PREFIX}MtGenomes_$LOCUS.fasta -o ref
			ls PrePaganAligned_*fasta | sed -e 's/Pre//' -e 's/\.fasta//' | parallel --no-notice -j $THREADS "pagan2 -a ref.fas -r ref.tre -q Pre{}.fasta -o {} --pileup --no-terminal-edges --silent"
			ls PaganAligned*fas | sed 's/\.fas//' | parallel --no-notice -j $THREADS "ntrlvdFasta2Fasta {}.fas {}.fasta FALSE"
			rm PaganAligned*fas
			ls PaganAligned*fasta | parallel --no-notice -j $THREADS "tail -n2 {}" | sed 's/\(.\)>/\1\n>/' > PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{15\}\(----------\)/\1NNNNNNNNNNNNNNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{14\}\(----------\)/\1NNNNNNNNNNNNNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{13\}\(----------\)/\1NNNNNNNNNNNNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{12\}\(----------\)/\1NNNNNNNNNNNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{11\}\(----------\)/\1NNNNNNNNNNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{10\}\(----------\)/\1NNNNNNNNNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{9\}\(----------\)/\1NNNNNNNNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{8\}\(----------\)/\1NNNNNNNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{7\}\(----------\)/\1NNNNNNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{6\}\(----------\)/\1NNNNNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{5\}\(----------\)/\1NNNNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{4\}\(----------\)/\1NNNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{3\}\(----------\)/\1NNN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{2\}\(----------\)/\1NN\2/g' PaganAligned_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{1\}\(----------\)/\1N\2/g' PaganAligned_RAD.fasta
			cat <(cat $(ls PaganAligned*fasta | head -n1) | head -n ${LinesInMtGenomes}) PaganAligned_RAD.fasta > ${PREFIX}ALL_masked_aligned_$LOCUS.fas

			mv ${PREFIX}ALL_masked_aligned_$LOCUS.fas ${PREFIX}ALL_masked_aligned_${LOCUS}.fasta
		fi
	echo ""; echo `date` Making sure all individuals have same number of nucleotides plus indels
		# count up nucleotides
		FILE=${PREFIX}ALL_masked_aligned_${LOCUS}.fasta
		maxNUC=$(cat $FILE | paste - - | awk -v col=2 'BEGIN{print "count", "lineNum"}{print gsub(/./,"",$col) "\t" NR}' | cut -f1 | sort -n | tail -1)
		echo "     There are $maxNUC nucleotides in the longest sequence"
		# save tsv file that quants the missing nucs
		paste <(cat $FILE | paste - - | awk -v col=2 'BEGIN{print "count", "lineNum"}{print gsub(/./,"",$col) "\t" NR}' | tail -n+2 | cut -f1) <(cat $FILE | paste - -) > ${FILE%.*}.tsv

		# make array of uniq numbers of missing deletions at the end of the sequences in the fasta file
		fixSEQ=($(awk -v x=$maxNUC '{ if ($1 < x) { print x-$1 }}' ${FILE%.*}.tsv | sort -n | uniq))
		# use sed to add - to ends of lines
		for i in ${fixSEQ[@]}; do
			echo "     $i indels being added to short sequences"
			NUC=$(($maxNUC - $i))
			insertSTRING=$(printf '%.0s-' $(seq 1 $i))
			#echo $NUC
			#echo $insertSTRING
			sed -i "s/^\(.\{$NUC\}\)$/\1$insertSTRING/" $FILE
		done
	#clean out sites with all gaps and 'n's, make sure to skip the lines with names
		sed -e '/>/!s/[nN]/\-/g' ${PREFIX}ALL_masked_aligned_$LOCUS.fasta | \
		seaview -convert -output_format fasta -o ${PREFIX}ALL_masked_aligned_clean_$LOCUS.fasta -del_gap_only_sites -
		seaview -convert -output_format nexus -o ${PREFIX}ALL_masked_aligned_clean_$LOCUS.nex ${PREFIX}ALL_masked_aligned_clean_$LOCUS.fasta
}
export -f alignLocusBySample

#function to make consensus sequences from aligned fasta files
mkMetaMitoGenomes(){
	#assign arguments to variables
		name=$1[@]
		local nontargetIDs=("${!name}")
		name=$2[@]
		local targetIDs=("${!name}")
		name=$3[@]
		local POPS=("${!name}")
		local PREFIX=$4
		local LOCUS=$5
		local THREADS=$6
		local cvgForCall=$7
		local nontargetNAME=$8
		local targetNAME=$9

	#split up fasta for making consensus seqs and name files according to seq id
		#split up fasta by sequence and make list of file names
			csplit --quiet --digits=4 --prefix=split${PREFIX}_masked_aligned_clean_$LOCUS.fasta ${PREFIX}ALL_masked_aligned_clean_$LOCUS.fasta "/^>/+0" "{*}"
			rm split${PREFIX}_masked_aligned_clean_$LOCUS.fasta0000
			local fileNAMES=($(ls split${PREFIX}_masked_aligned_clean_$LOCUS.fasta*))
			local seqNAMES=($(parallel -j $THREADS -k --no-notice "grep '^>' {} | sed 's/>//g' " ::: ${fileNAMES[@]}))
		#rename the files by sequence name
			parallel -j $THREADS -k --no-notice --link "mv {1} ${PREFIX}{2}_masked_aligned_clean_$LOCUS.fasta" ::: ${fileNAMES[@]} ::: ${seqNAMES[@]}
			parallel -j $THREADS -k --no-notice --link "echo {1} {2}" ::: ${fileNAMES[@]} ::: ${seqNAMES[@]}

	#make consensus sequence from outlier fish
		parallel -j $THREADS -k --no-notice "cat ${PREFIX}{}_masked_aligned_clean_$LOCUS.fasta" ::: ${nontargetIDs[@]} | sed 's/\(.\)>/\1\n>/' | grep -v '^cat:' | grep -v '^cat:' > ${PREFIX}NonTargetTaxon_masked_aligned_clean_$LOCUS.fasta
		Rscript consensusSeq.R ${PREFIX}NonTargetTaxon_masked_aligned_clean_$LOCUS.fasta ${PREFIX}NonTargetTaxon_masked_aligned_consensus_$LOCUS.fasta $cvgForCall ${nontargetNAME}_MetaMtGen

	#make consensus sequence from normal fish
		parallel -j $THREADS -k --no-notice "cat ${PREFIX}{}_masked_aligned_clean_$LOCUS.fasta" ::: ${targetIDs[@]} |  sed 's/\(.\)>/\1\n>/' | grep -v '^cat:' | grep -v '^cat:' > ${PREFIX}TargetTaxon_masked_aligned_clean_$LOCUS.fasta
		Rscript consensusSeq.R ${PREFIX}TargetTaxon_masked_aligned_clean_$LOCUS.fasta ${PREFIX}TargetTaxon_masked_aligned_consensus_$LOCUS.fasta $cvgForCall ${targetNAME}_MetaMtGen

	#make consensus sequence from normal fish for each sample location
		ncbiNAMES=($(echo ${seqNAMES[@]} | tr " " "\n" ))
		for i in ${POPS[@]}; do
			popIDs=($(echo ${targetIDs[@]} | tr " " "\n" | grep "^$i"))
			parallel -j $THREADS -k --no-notice "cat ${PREFIX}{}_masked_aligned_clean_$LOCUS.fasta" ::: ${popIDs[@]} | sed 's/\(.\)>/\1\n>/' | grep -v '^cat:' | grep -v '^cat:' > ${PREFIX}${i}_masked_aligned_clean_$LOCUS.fasta
			Rscript consensusSeq.R ${PREFIX}${i}_masked_aligned_clean_$LOCUS.fasta ${PREFIX}${i}_masked_aligned_consensus_$LOCUS.fasta $cvgForCall ${i}_MetaMtGen

			#get list of NCBI seqs
				ncbiNAMES=($(echo ${ncbiNAMES[@]} | tr " " "\n" | grep -v "$i"))
		done

	#get ncbi seqs
		parallel -j $THREADS -k --no-notice "cat ${PREFIX}{}_masked_aligned_clean_$LOCUS.fasta" ::: ${ncbiNAMES[@]} > ${PREFIX}NCBI_masked_aligned_clean_$LOCUS.fasta

	#concatenate consensus sequences into 1 file
		cat ${PREFIX}NCBI_masked_aligned_clean_$LOCUS.fasta ${PREFIX}*_masked_aligned_consensus_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_metamitogen_$LOCUS.fasta
	#mafft makes a bunch of sites with all indels due to a few poorly aligned sequences, remove gap only sites
		seaview -convert -output_format fasta -o ${PREFIX}ALL_masked_aligned_clean_metamitogen_$LOCUS.fasta -del_gap_only_sites ${PREFIX}ALL_masked_aligned_metamitogen_$LOCUS.fasta
	#convert fasta to nexus non interleaved
		seaview -convert -output_format nexus -o ${PREFIX}ALL_masked_aligned_clean_metamitogen_$LOCUS.nex ${PREFIX}ALL_masked_aligned_clean_metamitogen_$LOCUS.fasta
}

#function to maximize the number of bp retained at the expense of retaining individuals
maximizeBP(){
	#read arguments into variables
	local INFILE=$1
	local pctMissCall=$2
	Rscript maximizeBP.R ${INFILE%.*} $pctMissCall
	seaview -convert -output_format nexus -o ${INFILE%.*}_$pctMissCall.nex ${INFILE%.*}_$pctMissCall.fasta
}
