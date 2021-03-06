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

# function to determine if fasta file is interleaved
chkNTRLVD(){
	local INFILE=$1
	local numLines=$(wc -l $INFILE | tr -s ' ' '\t' | cut -f1)
	local numSeqs=$(grep -c '^>' $INFILE)
	local numLinesPerSeq=$(echo "$numLines / $numSeqs" | bc)
	if [ $(bc -l <<< "$numLinesPerSeq > 2") -eq 1 ]; then
		echo TRUE
	fi
}
export -f chkNTRLVD

chkINDELS(){
	local INFILE=$1
	local indelLines=$(grep -v '^>' $INFILE | grep -c '-')
	if [ $indelLines -gt 0 ]; then
		echo TRUE
	fi
}
export -f chkINDELS

chkSeqNAMES(){
	local INFILE=$1
	local numOffendingSeqNames=$(grep '^>' $INFILE | grep -c '[ |]')
	if [ $numOffendingSeqNames -gt 0 ]; then
		echo TRUE
	fi
}
export -f chkSeqNAMES


# function to convert an interleaved fasta to a normal fasta
ntrlvdFasta2Fasta(){
	if [ ! -z "$3" ] ; then
		local INFILE=$1
		local OUTFILE=$2
		local SORT=$3
		if [ "$SORT" == "TRUE" ]; then
			cat $INFILE | tr -d "\r" | sed 's/\(^>.*\)$/\1\t/' | tr -d "\n" | sed -e $'s/>/\\n>/g' | tail -n+2 | sort | tr "\t" "\n" > $OUTFILE
		else
			cat $INFILE | tr -d "\r" | sed 's/\(^>.*\)$/\1\t/' | tr -d "\n" | sed -e $'s/>/\\n>/g' | tail -n+2 | tr "\t" "\n" > $OUTFILE
		fi
	else
		local INFILE=/dev/stdin
		local SORT=$1
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
bam2GENO(){
		local ID=$1							#*_Pfa*
		local ID2=$1$2						#*_Pfa*$CUTOFFS-FLTR_F4
		local REF=$3						#unmasked reference genome
		
	# make output dir	
		local OUTDIR=out_bam2GENO

	# calculate coverage at each base from bam file and return a bed file with positions that have no coverage
		bedtools genomecov -ibam $ID2.bam -bga | grep -P '\t0$' > ./${OUTDIR}/$ID2.bed
		
	# create masked fasta for each individual
		bedtools maskfasta -fi $REF -fo ./${OUTDIR}/${ID2}_masked_ref.fasta -bed ./${OUTDIR}/$ID2.bed
		
	# make vcf files
		bcftools mpileup --threads 1 -d 30000 -q 30 -Q 20 -A -O z -o ./${OUTDIR}/${ID2}_masked_pile.vcf.gz -f ./${OUTDIR}/${ID2}_masked_ref.fasta ${ID2}.bam
		
	# call genotypes in vcf, force ploidy=haploid
		bcftools call --threads 1 -m --ploidy 1 -O z -o ./${OUTDIR}/${ID2}_masked_calls.vcf.gz ./${OUTDIR}/${ID2}_masked_pile.vcf.gz
		
	# combine snps and indels into 1 multiallelic call
		bcftools norm -f ./${OUTDIR}/${ID2}_masked_ref.fasta -m +any -O z -o ./${OUTDIR}/${ID2}_masked_calls_normalized.vcf.gz ./${OUTDIR}/${ID2}_masked_calls.vcf.gz
		tabix ./${OUTDIR}/${ID2}_masked_calls_normalized.vcf.gz
		
	# generate 1 consensus for each fishSample
		bcftools consensus -s $ID -M 'N' -f ./${OUTDIR}/${ID2}_masked_ref.fasta -o ./${OUTDIR}/${ID2}_masked_consensus.fasta ./${OUTDIR}/${ID2}_masked_calls_normalized.vcf.gz 
		
	# add ID to name of consensus
		sed -i "s/^>/>${ID}_/g" ./${OUTDIR}/${ID2}_masked_consensus.fasta	
}
export -f bam2GENO

#function to get locus from masked consensus sequences, mito genomes, and NCBI nucleotide records, clean and align
aliGENO(){

	# assign arguments to variables
		local PREFIX=$1
		local THREADS=$2
		local midFILE=$3
		local POSITIONS=$4
		local LOCUS=$5
		local mtGenPATTERN="$6"	
		local LONGALIGNMENT=$7
		local REF=$8
		local GENBANK=$9


		# local name=$7[@]
		# local IDs=("${!name}")
	
	echo ""; echo PREFIX: $PREFIX
	echo THREADS: $THREADS
	echo midFILE: $midFILE
	echo POSITIONS: $POSITIONS
	echo LOCUS: $LOCUS
	echo mtGenPATTERN: $mtGenPATTERN
	echo Long alignment: $LONGALIGNMENT
	echo Genbank: $GENBANK
	echo REF: $REF; echo ""
	
	# set input and make output dir
	if [ -d out_bam2GENO ]; then
		INDIR=out_bam2GENO
	else
		INDIR=.
	fi
	OUTDIR=out_aliGENO
	if [ -d ${OUTDIR} ]; then
		echo ""; echo `date` $OUTDIR EXISTS, REMOVING CONTENTS FROM PREVIOUS RUN...
		rm -rf ${OUTDIR}
	fi
	mkdir ${OUTDIR}

	#get locus from all individuals and align, $POSITIONS determines the locus
	echo ""; echo `date` EXTRACTING POSITIONS $POSITIONS FOR ALIGNMENT...
		echo ${IDs[@]} | tr " " "\n" | parallel -j $THREADS -k --no-notice "echo \>{} && tail -n +2 ${INDIR}/{}${midFILE}_masked_consensus.fasta | tr '\n' '\t' | sed 's/\t//g' | sed 's/ */\t/g' | cut -f $POSITIONS | sed 's/\t//g' " > ./${OUTDIR}/${PREFIX}RAD_masked_$LOCUS.fasta
		
	#remove individuals with no sequences or all N
	echo ""; echo `date` REMOVING INDIVIDUALS WITH NO NUCLEOTIDES CALLED...
		sed 's/^NN*$//g' ./${OUTDIR}/${PREFIX}RAD_masked_$LOCUS.fasta | grep -P -B 1 '^[ACTGN@]+$' | grep -v '\-\-' > ./${OUTDIR}/${PREFIX}RAD_masked_${LOCUS}_clean.fasta
	
	#make fasta from mtGenome sequences
	if [ ! -z "$mtGenPATTERN" ]; then
		echo ""; echo `date` GATHERING ALL mtGENOMEs WITH PATTERN=$mtGenPATTERN AND EXTRACTING POSITIONS $POSITIONS
		ls $mtGenPATTERN | sed 's/.fasta//g' | parallel -j $THREADS -k --no-notice "echo \>{} && tail -n +2 {}.fasta | tr '\n' '\t' | sed 's/\t//g' | sed 's/ */\t/g' | cut -f $POSITIONS | sed 's/\t//g' " > ./${OUTDIR}/${PREFIX}MtGenomes_$LOCUS.fasta
		LinesInMtGenomes=$(wc -l ./${OUTDIR}/${PREFIX}MtGenomes_$LOCUS.fasta | cut -d" " -f1)
		NumMtGenomes=$(grep -c "^>" ./${OUTDIR}/${PREFIX}MtGenomes_$LOCUS.fasta)
	else
		echo ""; echo `date` NO mtGENOMES SPECIFIED, CONTINUING WITHOUT THEM...
	fi
	#Make fasta from other NCBI nucleotide sequences
		#need to code, right now, just download as 1 big fasta

	#Combine fastas, make sure to add genbank nucleotide recs as neccessary
	# Fix GENBANKFASTA file
	if [ ! -z "$GENBANK" ]; then
		echo ""; echo `date` ADDING $GENBANK SEQUENCES FROM GENBANK...
		#remove carriage returns
		sed -i 's/\r//g' $GENBANK
		if [ "$(chkNTRLVD $GENBANK)" == "TRUE" ]; then
			echo ""; echo `date` $GENBANK IS INTERLEAVED FASTA, CREATING SEQUENTIAL FASTA...
				local sqntlGENBANK=${GENBANK%.*}_sqntl.${GENBANK##*.}
				ntrlvdFasta2Fasta $GENBANK $sqntlGENBANK FALSE
				GENBANK=$sqntlGENBANK
		fi
		if [ "$(chkINDELS $GENBANK)" == "TRUE" ]; then
			echo ""; echo `date` $GENBANK HAS INDELS, REMOVING INDELS...
				local noindelGENBANK=${GENBANK%.*}_noIndel.${GENBANK##*.}
				paste <(cat $GENBANK | paste - - | cut -f1) \
					<(cat $GENBANK | paste - - | cut -f2 | sed 's/-//g') | \
					tr '\t' '\n' > $noindelGENBANK
				GENBANK=$noindelGENBANK
		fi
		if [ "$(chkSeqNAMES $GENBANK)" == "TRUE" ]; then
			echo ""; echo      `date` $GENBANK HAS SEQUENCE NAMES WITH \| AND SPACES, RENAMING SEQUENCES...
				local cleanedGENBANK=${GENBANK%.*}_cleaned.${GENBANK##*.}
				paste <(cat $GENBANK | paste - - | cut -f1 | sed -e 's/|.*$//' -e 's/ .*$//') \
					<(cat $GENBANK | paste - - | cut -f2) | \
					tr '\t' '\n' > $cleanedGENBANK
				GENBANK=$cleanedGENBANK
		fi
	else
		echo ""; echo `date` NO ADDITIONAL GENBANK SEQUENCES SPECIFIED, CONTINUING WITHOUT THEM...
	fi
	
	#echo ""; echo `date` CONCATENATING FASTAs...
	
	#cat ./${OUTDIR}/${PREFIX}MtGenomes_$LOCUS.fasta > ./${OUTDIR}/test1.fasta
	#cat ./${OUTDIR}/${PREFIX}RAD_masked_${LOCUS}_clean.fasta > ./${OUTDIR}/test2.fasta
	#cat $GENBANK  > ./${OUTDIR}/test3.fasta

	# Align all_*.fasta
	echo ""; echo `date` ALIGNING SEQUENCES WITH $(if [ "$LONGALIGNMENT" == "TRUE" ]; then echo -n MAFFT; else echo -n pagan2; fi)...
		echo ""
		if [ "$LONGALIGNMENT" == "TRUE" ] || [ "$LONGALIGNMENT" == "T" ]; then
			cat ./${OUTDIR}/${PREFIX}MtGenomes_$LOCUS.fasta ./${OUTDIR}/${PREFIX}RAD_masked_${LOCUS}_clean.fasta $GENBANK  > ./${OUTDIR}/${PREFIX}ALL_masked_$LOCUS.fasta
			#clustalw -infile=${PREFIX}ALL_masked_$LOCUS.fasta -align -type=DNA -output=NEXUS -outfile=${PREFIX}ALL_masked_aligned_$LOCUS.nex
			#clustalo -infile=${PREFIX}ALL_masked_$LOCUS.fasta -t DNA --outfmt=fa -outfile=${PREFIX}ALL_masked_aligned_$LOCUS.fasta --threads $THREADS
			#mafft --thread $THREADS ${PREFIX}ALL_masked_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_$LOCUS.fasta
			#mafft --thread $THREADS --ep 0.123 ${PREFIX}ALL_masked_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_$LOCUS.fasta
			mafft --thread $THREADS --globalpair --maxiterate 1000 ./${OUTDIR}/${PREFIX}ALL_masked_$LOCUS.fasta > ./${OUTDIR}/${PREFIX}ALL_masked_aligned_$LOCUS.fasta
		else
			# my solution for pagan alignment: align 1 individual at a time to the mulitiple ref genomes
			echo `date` ALIGNING RAD DATA TO MITOGENOMES...
			cat ./${OUTDIR}/${PREFIX}RAD_masked_${LOCUS}_clean.fasta $GENBANK  > ./${OUTDIR}/${PREFIX}ALL_masked_$LOCUS.fasta
			sed 's/N\{20,\}/N/g' ./${OUTDIR}/${PREFIX}ALL_masked_$LOCUS.fasta > ./${OUTDIR}/${PREFIX}ALL_masked_${LOCUS}_clean.fasta

			# causes error in the next parallel command if there are too many samples
			#IndivNames=($(cat ./${OUTDIR}/${PREFIX}ALL_masked_${LOCUS}_clean.fasta | paste - - | sed 's/^>//' | cut -f1))
			#echo ""; echo Print IndivNames:; echo ${IndivNames[@]}; echo ""
			#IndivSeqs=($(cat ./${OUTDIR}/${PREFIX}ALL_masked_${LOCUS}_clean.fasta | paste - - | sed 's/^>//' | cut -f2))

			#solution to error
			cat ./${OUTDIR}/${PREFIX}ALL_masked_${LOCUS}_clean.fasta | paste - - | sed 's/^>//' | cut -f1 > ./${OUTDIR}/IndivNames.txt
			cat ./${OUTDIR}/${PREFIX}ALL_masked_${LOCUS}_clean.fasta | paste - - | sed 's/^>//' | cut -f2 > ./${OUTDIR}/IndivSeqs.txt

			# parallel --no-notice --link -j $THREADS "printf '>{1}\n{2}\n' > ./${OUTDIR}/PrePaganAligned_{1}.fasta " ::: ${IndivNames[@]} ::: ${IndivSeqs[@]}
			parallel --no-notice --link -j $THREADS "printf '>{1}\n{2}\n' > ./${OUTDIR}/PrePaganAligned_{1}.fasta " :::: ./${OUTDIR}/IndivNames.txt :::: ./${OUTDIR}/IndivSeqs.txt

			# error with 1 ref genome
			pagan2 -s ./${OUTDIR}/${PREFIX}MtGenomes_$LOCUS.fasta -o ./${OUTDIR}/ref

			if [ $NumMtGenomes -eq 1 ]; then
				cat ./${OUTDIR}/${PREFIX}MtGenomes_$LOCUS.fasta > ./${OUTDIR}/ref.fas
			fi

			ls ./${OUTDIR}/PrePaganAligned_*fasta | sed -e "s/\.\/${OUTDIR}\///" -e 's/Pre//' -e 's/\.fasta//' | parallel --no-notice -j $THREADS "pagan2 -a ./${OUTDIR}/ref.fas -r ./${OUTDIR}/ref.tre -q ./${OUTDIR}/Pre{}.fasta -o ./${OUTDIR}/{} --pileup-alignment --no-terminal-edges --silent"

			ls ./${OUTDIR}/PaganAligned*fas | sed -e "s/\.\/${OUTDIR}\///" -e 's/\.fas//' | parallel --no-notice -j $THREADS "ntrlvdFasta2Fasta ./${OUTDIR}/{}.fas ./${OUTDIR}/{}.fasta FALSE"

			ls ./${OUTDIR}/PaganAligned*fasta | parallel --no-notice -j $THREADS "tail -n2 {}" | sed 's/\(.\)>/\1\n>/' > ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{15\}\(----------\)/\1NNNNNNNNNNNNNNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{14\}\(----------\)/\1NNNNNNNNNNNNNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{13\}\(----------\)/\1NNNNNNNNNNNNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{12\}\(----------\)/\1NNNNNNNNNNNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{11\}\(----------\)/\1NNNNNNNNNNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{10\}\(----------\)/\1NNNNNNNNNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{9\}\(----------\)/\1NNNNNNNNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{8\}\(----------\)/\1NNNNNNNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{7\}\(----------\)/\1NNNNNNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{6\}\(----------\)/\1NNNNNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{5\}\(----------\)/\1NNNNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{4\}\(----------\)/\1NNNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{3\}\(----------\)/\1NNN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{2\}\(----------\)/\1NN\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			sed -i 's/\(----------\)[ACTGN]\{1\}\(----------\)/\1N\2/g' ./${OUTDIR}/PaganAlign_RAD.fasta
			cat <(cat $(ls ./${OUTDIR}/PaganAligned*fasta | head -n1) | head -n ${LinesInMtGenomes}) ./${OUTDIR}/PaganAlign_RAD.fasta > ./${OUTDIR}/${PREFIX}ALL_masked_aligned_$LOCUS.fasta
		fi
		
	echo ""; echo `date` Making sure all individuals have same number of nucleotides plus indels
		# count up nucleotides
		FILE=./${OUTDIR}/${PREFIX}ALL_masked_aligned_${LOCUS}.fasta
		maxNUC=$(cat $FILE | paste - - | awk -v col=2 'BEGIN{print "count", "lineNum"}{print gsub(/./,"",$col) "\t" NR}' | cut -f1 | sort -n | tail -1)
		echo "     There are $maxNUC nucleotides in the longest sequence"
		# save tsv file that quants the missing nucs
		paste <(cat $FILE | paste - - | awk -v col=2 'BEGIN{print "count", "lineNum"}{print gsub(/./,"",$col) "\t" NR}' | tail -n+2 | cut -f1) <(cat $FILE | paste - -) > ${FILE%.*}.tsv

		# make array of uniq numbers of missing deletions at the end of the sequences in the fasta file
		# this is slow
		#fixSEQ=($(awk -v x=$maxNUC '{ if ($1 < x) { print x-$1 }}' ${FILE%.*}.tsv | sort -n | uniq))
		# use sed to add - to ends of lines
		#for i in ${fixSEQ[@]}; do
		#	echo "     $i indels being added to short sequences"
		#	NUC=$(($maxNUC - $i))
		#	insertSTRING=$(printf '%.0s-' $(seq 1 $i))
		#	#echo $NUC
		#	#echo $insertSTRING
		#	sed -i "s/^\(.\{$NUC\}\)$/\1$insertSTRING/" $FILE
		#done
		# alternative sed that should be faster
		fixSEQ=($(awk -v x=$maxNUC '{ if ($1 < x) { print x-$1 } else { print 0 } }' ${FILE%.*}.tsv ))
		cp ${FILE%.*}.tsv ${FILE%.*}_eqLen.tsv
		for i in ${!fixSEQ[@]}; do
			# echo "     ${fixSEQ[$i]} indels being added to short sequences"
			if [ ${fixSEQ[$i]} == 0 ]; then
				insertSTRING=""
			else
				insertSTRING=$(printf '%.0s-' $(seq 1 ${fixSEQ[$i]}))
			fi
				#echo $NUC
				#echo $insertSTRING
			j=$(($i + 1))
			sed -i "$j s/$/$insertSTRING/" ${FILE%.*}_eqLen.tsv
		done
		cat ${FILE%.*}_eqLen.tsv | cut -f 2-3 | tr "\t" "\n" > ${FILE%.*}_eqLen.fasta

	# clean out sites with all gaps and 'n's, make sure to skip the lines with names
		#sed -e '/>/!s/[nN]/\-/g' ./${OUTDIR}/${PREFIX}ALL_masked_aligned_$LOCUS.fasta | \
		sed -e '/>/!s/[nN]/\-/g' ${FILE%.*}_eqLen.fasta | \
			seaview -convert -output_format fasta -o ./${OUTDIR}/${PREFIX}ALL_masked_aligned_clean_$LOCUS.fasta -del_gap_only_sites -
		seaview -convert -output_format nexus -o ./${OUTDIR}/${PREFIX}ALL_masked_aligned_clean_$LOCUS.nex ./${OUTDIR}/${PREFIX}ALL_masked_aligned_clean_$LOCUS.fasta
}		
export -f aliGENO

#function to make consensus sequences from aligned fasta files
mkMETAGENO(){
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
		
	# make output dir
	INDIR=out_aliGENO
	if [ -f out_aliGENO ]; then
		INDIR=out_aliGENO
	else
		INDIR=.
	fi

	OUTDIR=out_metaGENO
	if [ -d ${OUTDIR} ]; then
		rm -rf ${OUTDIR}
	fi
	mkdir ${OUTDIR}
	
	#split up fasta for making consensus seqs and name files according to seq id
		#split up fasta by sequence and make list of file names
			csplit --quiet --digits=4 --prefix=split${PREFIX}_masked_aligned_clean_$LOCUS.fasta ${INDIR}/${PREFIX}ALL_masked_aligned_clean_$LOCUS.fasta "/^>/+0" "{*}"
			rm split${PREFIX}_masked_aligned_clean_$LOCUS.fasta0000
			mv split${PREFIX}_masked_aligned_clean_$LOCUS.fasta* ${OUTDIR}
			local fileNAMES=($(ls ./${OUTDIR}/split${PREFIX}_masked_aligned_clean_$LOCUS.fasta*))
			local seqNAMES=($(parallel -j $THREADS -k --no-notice "grep '^>' {} | sed 's/>//g' " ::: ${fileNAMES[@]}))
		#rename the files by sequence name
			parallel -j $THREADS -k --no-notice --link "mv {1} ./${OUTDIR}/${PREFIX}{2}_masked_aligned_clean_$LOCUS.fasta" ::: ${fileNAMES[@]} ::: ${seqNAMES[@]}
			# parallel -j $THREADS -k --no-notice --link "echo {1} {2}" ::: ${fileNAMES[@]} ::: ${seqNAMES[@]}

	#make consensus sequence from outlier fish
		parallel -j $THREADS -k --no-notice "cat ./${OUTDIR}/${PREFIX}{}_masked_aligned_clean_$LOCUS.fasta 2> /dev/null" ::: ${nontargetIDs[@]} | sed 's/\(.\)>/\1\n>/' | grep -v '^cat:' | grep -v '^cat:' > ./${OUTDIR}/${PREFIX}NonTargetTaxon_masked_aligned_clean_$LOCUS.fasta
		consensusSEQ.R ./${OUTDIR}/${PREFIX}NonTargetTaxon_masked_aligned_clean_$LOCUS.fasta ./${OUTDIR}/${PREFIX}NonTargetTaxon_masked_aligned_metageno_$LOCUS.fasta $cvgForCall ${nontargetNAME}_MetaMtGen

	#make consensus sequence from normal fish
		parallel -j $THREADS -k --no-notice "cat ./${OUTDIR}/${PREFIX}{}_masked_aligned_clean_$LOCUS.fasta 2> /dev/null" ::: ${targetIDs[@]} |  sed 's/\(.\)>/\1\n>/' | grep -v '^cat:' | grep -v '^cat:' > ./${OUTDIR}/${PREFIX}TargetTaxon_masked_aligned_clean_$LOCUS.fasta
		consensusSEQ.R ./${OUTDIR}/${PREFIX}TargetTaxon_masked_aligned_clean_$LOCUS.fasta ./${OUTDIR}/${PREFIX}TargetTaxon_masked_aligned_metageno_$LOCUS.fasta $cvgForCall ${targetNAME}_MetaMtGen

	#make consensus sequence from normal fish for each sample location
		ncbiNAMES=($(echo ${seqNAMES[@]} | tr " " "\n" ))
		for i in ${POPS[@]}; do
			popIDs=($(echo ${targetIDs[@]} | tr " " "\n" | grep "^$i"))
			parallel -j $THREADS -k --no-notice "cat ./${OUTDIR}/${PREFIX}{}_masked_aligned_clean_$LOCUS.fasta 2> /dev/null" ::: ${popIDs[@]} | sed 's/\(.\)>/\1\n>/' | grep -v '^cat:' | grep -v '^cat:' > ./${OUTDIR}/${PREFIX}${i}-POP_masked_aligned_clean_$LOCUS.fasta
			consensusSEQ.R ./${OUTDIR}/${PREFIX}${i}-POP_masked_aligned_clean_$LOCUS.fasta ./${OUTDIR}/${PREFIX}${i}-POP_masked_aligned_metageno_$LOCUS.fasta $cvgForCall ${i}_MetaMtGen

			#get list of NCBI seqs
				ncbiNAMES=($(echo ${ncbiNAMES[@]} | tr " " "\n" | grep -v "$i"))
		done

	#get ncbi seqs
		parallel -j $THREADS -k --no-notice "cat ./${OUTDIR}/${PREFIX}{}_masked_aligned_clean_$LOCUS.fasta 2> /dev/null" ::: ${ncbiNAMES[@]} > ./${OUTDIR}/${PREFIX}NCBI_masked_aligned_clean_$LOCUS.fasta

	#concatenate consensus sequences into 1 file
		cat ./${OUTDIR}/${PREFIX}NCBI_masked_aligned_clean_$LOCUS.fasta ./${OUTDIR}/${PREFIX}*_masked_aligned_metageno_$LOCUS.fasta > ./${OUTDIR}/${PREFIX}ALL_masked_aligned_metageno_$LOCUS.fasta
		cat ./${OUTDIR}/${PREFIX}NCBI_masked_aligned_clean_$LOCUS.fasta ./${OUTDIR}/${PREFIX}*Taxon_masked_aligned_metageno_$LOCUS.fasta > ./${OUTDIR}/${PREFIX}TAXA_masked_aligned_metageno_$LOCUS.fasta
		cat ./${OUTDIR}/${PREFIX}NCBI_masked_aligned_clean_$LOCUS.fasta ./${OUTDIR}/${PREFIX}*POP_masked_aligned_metageno_$LOCUS.fasta > ./${OUTDIR}/${PREFIX}POPS_masked_aligned_metageno_$LOCUS.fasta
	#mafft makes a bunch of sites with all indels due to a few poorly aligned sequences, remove gap only sites
		seaview -convert -output_format fasta -o ./${OUTDIR}/${PREFIX}ALL_masked_aligned_clean_metageno_$LOCUS.fasta -del_gap_only_sites ./${OUTDIR}/${PREFIX}ALL_masked_aligned_metageno_$LOCUS.fasta
		seaview -convert -output_format fasta -o ./${OUTDIR}/${PREFIX}TAXA_masked_aligned_clean_metageno_$LOCUS.fasta -del_gap_only_sites ./${OUTDIR}/${PREFIX}TAXA_masked_aligned_metageno_$LOCUS.fasta
		seaview -convert -output_format fasta -o ./${OUTDIR}/${PREFIX}POPS_masked_aligned_clean_metageno_$LOCUS.fasta -del_gap_only_sites ./${OUTDIR}/${PREFIX}POPS_masked_aligned_metageno_$LOCUS.fasta
#convert fasta to nexus non interleaved
		seaview -convert -output_format nexus -o ./${OUTDIR}/${PREFIX}ALL_masked_aligned_clean_metageno_$LOCUS.nex ./${OUTDIR}/${PREFIX}ALL_masked_aligned_clean_metageno_$LOCUS.fasta
		seaview -convert -output_format nexus -o ./${OUTDIR}/${PREFIX}TAXA_masked_aligned_clean_metageno_$LOCUS.nex ./${OUTDIR}/${PREFIX}TAXA_masked_aligned_clean_metageno_$LOCUS.fasta
		seaview -convert -output_format nexus -o ./${OUTDIR}/${PREFIX}POPS_masked_aligned_clean_metageno_$LOCUS.nex ./${OUTDIR}/${PREFIX}POPS_masked_aligned_clean_metageno_$LOCUS.fasta
}

#function to maximize the number of bp retained at the expense of retaining individuals
fltrGENOSITES(){
	#read arguments into variables
	local INFILE=$1
	local pctMissCall=$2
	fltrGENOSITES.R ${INFILE%.*} $pctMissCall
	seaview -convert -output_format nexus -o ${INFILE%.*}_$pctMissCall.nex ${INFILE%.*}_$pctMissCall.fasta
}
