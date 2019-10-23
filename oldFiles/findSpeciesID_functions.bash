#!/bin/bash

#this script 

#user defined variables
THREADS=7
CUTOFFS=".48.48"							#dDocent cutoffs used for reference genome
FILTER=4									#filter reads with this samtools flag
PREFIX=fish									#prefix on files created
POSITIONS=5464-6500							#start and end positions of mtDNA fragment to excise
LOCUS=begCOI								#name of locus
mtGENOMEs="*Puntioplites*mtGenome.fasta"	#ls search string to return mitoGenomes, names should end with _mtGenome.fasta
OUTLIERS=RAD_OUTLIER_Pfalcifer_fish.txt		#name of file containing fish in outlier group

#automatic variables
REF=reference${CUTOFFS}.fasta
BAMLIST=bamlist${CUTOFFS}.list
catFILE1=cat${CUTOFFS}-FLTR_F4.bam 
bamLIST2=bamlist${CUTOFFS}.list2
catFILE2=cat${CUTOFFS}-FLTR_F$FILTER.bam 


#function to create consensus sequences using individually masked reference genomes
	#superior to using a single unmasked reference genome
	#designed to be run in parallel
#mkMaskedCONSENSI(){
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
export -f mkMaskedCONSENSI

#function to get locus from masked consensus sequences, mito genomes, and NCBI nucleotide records, clean and align
alignLocusBySample(){
	#assign arguments to variables
		local PREFIX=$1
		local THREADS=$2
		local midFILE=$3
		local POSITIONS=$4
		local LOCUS=$5
		#eval mtgenomes="$6"		#cant figure this one out
		local GENBANK=$6
		name=$7[@]
		local IDs=("${!name}")

		
	#get locus from all individuals and align, $POSITIONS determines the locus
		echo ${IDs[@]} | tr " " "\n" | parallel -j $THREADS -k --no-notice "echo \>{} && tail -n +2 {}${midFILE}_masked_consensus.fasta | tr '\n' '\t' | sed 's/\t//g' | sed 's/ */\t/g' | cut -f $POSITIONS | sed 's/\t//g' " > ${PREFIX}RAD_masked_$LOCUS.fasta
	#remove individuals with no sequences or all N
		sed 's/^NN*$//g' ${PREFIX}RAD_masked_$LOCUS.fasta | grep -P -B 1 '^[ACTGN@]+$' | grep -v '\-\-' > ${PREFIX}RAD_masked_${LOCUS}_clean.fasta
	#make fasta from mtGenome sequences
		ls $mtGENOMEs | sed 's/_mtGenome.fasta//g' | parallel -j $THREADS -k --no-notice "echo \>{} && tail -n +2 {}_mtGenome.fasta | tr '\n' '\t' | sed 's/\t//g' | sed 's/ */\t/g' | cut -f $POSITIONS | sed 's/\t//g' " > ${PREFIX}MtGenomes_$LOCUS.fasta
	#Make fasta from other NCBI nucleotide sequences
		#need to code, right now, just download as 1 big fasta
	#Combine fastas, make sure to add genbank nucleotide recs as neccessary
		echo $GENBANK
		cat ${PREFIX}MtGenomes_$LOCUS.fasta ${PREFIX}RAD_masked_${LOCUS}_clean.fasta $GENBANK  > ${PREFIX}ALL_masked_$LOCUS.fasta 
	#Align all_COI.fasta
		#clustalw -infile=${PREFIX}ALL_masked_$LOCUS.fasta -align -type=DNA -output=NEXUS -outfile=${PREFIX}ALL_masked_aligned_$LOCUS.nex
		#clustalo -infile=${PREFIX}ALL_masked_$LOCUS.fasta -t DNA --outfmt=fa -outfile=${PREFIX}ALL_masked_aligned_$LOCUS.fasta --threads $THREADS 
		#mafft --thread $THREADS ${PREFIX}ALL_masked_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_$LOCUS.fasta
		#mafft --thread $THREADS --ep 0.123 ${PREFIX}ALL_masked_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_$LOCUS.fasta
		mafft --thread $THREADS --globalpair --maxiterate 1000 ${PREFIX}ALL_masked_$LOCUS.fasta > ${PREFIX}ALL_masked_aligned_$LOCUS.fasta
	#clean out sites with all gaps and 'n's, make sure to skip the lines with names
		sed -e '/>/!s/n/\-/g' ${PREFIX}ALL_masked_aligned_$LOCUS.fasta | \
		#seaview -convert -output_format nexus -o ${PREFIX}ALL_masked_aligned_$LOCUS.nex ${PREFIX}ALL_masked_aligned_$LOCUS.fasta 
		seaview -convert -output_format nexus -o ${PREFIX}ALL_masked_aligned_clean_$LOCUS.nex -del_gap_only_sites -

}

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
		# #replace indels with U and N with indels so that consensus will return nucleotide if 2 or more present
			# sed -i '/>/!s/-/U/g' ${PREFIX}OUTLIERS_masked_aligned_$LOCUS.fasta
			# sed -i '/>/!s/N/-/g' ${PREFIX}OUTLIERS_masked_aligned_$LOCUS.fasta
		#cons --sequence ${PREFIX}OUTLIERS_masked_aligned_$LOCUS.fasta -outseq ${PREFIX}OUTLIERS_masked_aligned_consensus_$LOCUS.fasta -plurality 0.5 -identity 2 -name Olier_Consens
		Rscript consensusSeq.R ${PREFIX}OUTLIERS_masked_aligned_$LOCUS.fasta ${PREFIX}OUTLIERS_masked_aligned_consensus_$LOCUS.fasta $cvgForCall Outlier_Consensus
		
		# #replace U with indels and change back the rest
			# sed -i '/>/!s/-/N/g' ${PREFIX}OUTLIERS_masked_aligned_consensus_$LOCUS.fasta
			# sed -i '/>/!s/[Uu]/-/g' ${PREFIX}OUTLIERS_masked_aligned_consensus_$LOCUS.fasta
			# sed -i '/>/!s/-/N/g' ${PREFIX}OUTLIERS_masked_aligned_$LOCUS.fasta
			# sed -i '/>/!s/U/-/g' ${PREFIX}OUTLIERS_masked_aligned_$LOCUS.fasta
			
	#make consensus sequence from normal fish
		parallel -j $THREADS -k --no-notice "cat ${PREFIX}{}_masked_aligned_$LOCUS.fasta" ::: ${normIDs[@]} | grep -v '^cat:' | grep -v '^cat:' > ${PREFIX}NORMALS_masked_aligned_$LOCUS.fasta
		# #replace indels with U and N with indels so that consensus will return nucleotide if 2 or more present
			# sed -i '/>/!s/-/U/g' ${PREFIX}NORMALS_masked_aligned_$LOCUS.fasta
			# sed -i '/>/!s/N/-/g' ${PREFIX}NORMALS_masked_aligned_$LOCUS.fasta
			#cons --sequence ${PREFIX}NORMALS_masked_aligned_$LOCUS.fasta -outseq ${PREFIX}NORMALS_masked_aligned_consensus_$LOCUS.fasta -plurality 0.5 -identity 2 -name Normal_Consens
		Rscript consensusSeq.R ${PREFIX}NORMALS_masked_aligned_$LOCUS.fasta ${PREFIX}NORMALS_masked_aligned_consensus_$LOCUS.fasta $cvgForCall Normal_Consensus

		# #replace U with indels and change back the rest
			# sed -i '/>/!s/-/N/g' ${PREFIX}NORMALS_masked_aligned_consensus_$LOCUS.fasta
			# sed -i '/>/!s/[Uu]/-/g' ${PREFIX}NORMALS_masked_aligned_consensus_$LOCUS.fasta
			# sed -i '/>/!s/-/N/g' ${PREFIX}NORMALS_masked_aligned_$LOCUS.fasta
			# sed -i '/>/!s/U/-/g' ${PREFIX}NORMALS_masked_aligned_$LOCUS.fasta
			

	#make consensus sequence from normal fish for each sample location
		ncbiNAMES=($(echo ${seqNAMES[@]} | tr " " "\n" ))
		for i in ${SITES[@]}; do
			siteIDs=($(echo ${normIDs[@]} | tr " " "\n" | grep "^$i"))
			parallel -j $THREADS -k --no-notice "cat ${PREFIX}{}_masked_aligned_$LOCUS.fasta" ::: ${siteIDs[@]} | grep -v '^cat:' | grep -v '^cat:' > ${PREFIX}${i}_masked_aligned_$LOCUS.fasta
			# #replace indels with U and N with indels so that consensus will return nucleotide if 2 or more present
				# sed -i '/>/!s/-/U/g' ${PREFIX}${i}_masked_aligned_$LOCUS.fasta
				# sed -i '/>/!s/N/-/g' ${PREFIX}${i}_masked_aligned_$LOCUS.fasta
			#cons --sequence ${PREFIX}${i}_masked_aligned_$LOCUS.fasta -outseq ${PREFIX}${i}_masked_aligned_consensus_$LOCUS.fasta -plurality 0.5 -identity 2 -name ${i}_Consens
			Rscript consensusSeq.R ${PREFIX}${i}_masked_aligned_$LOCUS.fasta ${PREFIX}${i}_masked_aligned_consensus_$LOCUS.fasta $cvgForCall ${i}_Consensus

			# #replace U with indels and change back the rest
				# sed -i '/>/!s/-/N/g' ${PREFIX}${i}_masked_aligned_consensus_$LOCUS.fasta
				# sed -i '/>/!s/[Uu]/-/g' ${PREFIX}${i}_masked_aligned_consensus_$LOCUS.fasta
				# sed -i '/>/!s/-/N/g' ${PREFIX}${i}_masked_aligned_$LOCUS.fasta
				# sed -i '/>/!s/U/-/g' ${PREFIX}${i}_masked_aligned_$LOCUS.fasta
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
		seqret -sequence ${PREFIX}ALL_masked_aligned_clean_consensi_$LOCUS.fasta -outseq nexusnon::${PREFIX}ALL_masked_aligned_clean_consensi_$LOCUS.nex
}







#function to merge all bam files
catBAM(){
	#assign arguments to variables
		local CUTOFFS=$1
		local BAMLIST=bamlist${1}.list
		local catFILE1=cat${1}-FLTR_F4.bam
	#combine bam files from all individuals and index the concatenated bam file
		ls *_*${CUTOFFS}-FLTR_F4.bam > $BAMLIST
		samtools merge -@$THREADS -b $BAMLIST -f $catFILE1 &>/dev/null
		samtools index $catFILE1
	#evaluate mapping stats from bam files
		samtools flagstat $catFILE1 > cat${CUTOFFS}-FLTR_F4.flagstat
}


#function to merge bam files by group, uses RAD_OUTLIER_Pfalcifer_fish.txt
grpCatBAM(){
	#assign arguments to variables
		local CUTOFFS=$1
		local OUTLIERS=$2
		local BAMLIST=bamlist${1}.list
		local OUTLIST=outlier_$BAMLIST
		local NORMLIST=normal_$BAMLIST
		local outCatFILE=outCat${1}-FLTR_F4.bam
		local normCatFILE=normCat${1}-FLTR_F4.bam
	#combine bam files from all individuals and index the concatenated bam file
		ls *_*${CUTOFFS}-FLTR_F4.bam > $BAMLIST
		cat $OUTLIERS | parallel -j $THREADS -k --no-notice "grep {} $BAMLIST" > $OUTLIST
		grep -v -f $OUTLIERS $BAMLIST > $NORMLIST
		samtools merge -@$THREADS -b $NORMLIST -f $normCatFILE 
		samtools merge -@$THREADS -b $OUTLIST -f $outCatFILE 
		samtools index $outCatFILE
		samtools index $normCatFILE
	#evaluate mapping stats from bam files
		samtools flagstat $outCatFILE > outCat${CUTOFFS}-FLTR_F4.flagstat
		samtools flagstat $normCatFILE > normCat${CUTOFFS}-FLTR_F4.flagstat
}


#function to make consensus mtGenome for each sample
getConsensusBySample(){
	#assign arguments to variables
		local ref=reference${1}.fasta
		local catfile2=cat${1}-FLTR_F$2.bam
		local prefix=$3
		
	#make vcf file
		bcftools mpileup --threads $THREADS -d 20000 -q 30 -Q 20 -A -O z -o ${prefix}Pile.vcf.gz -f $ref $catfile2
	#call genotypes in vcf, force ploidy=haploid
		bcftools call --threads $THREADS -m --ploidy 1 -O z -o ${prefix}Calls.vcf.gz ${prefix}Pile.vcf.gz
	#combine snps and indels into 1 multiallelic call and index
		bcftools norm -f $ref -m +any -O z -o ${prefix}Calls_normalized.vcf.gz ${prefix}Calls.vcf.gz
		tabix ${prefix}Calls_normalized.vcf.gz
	#compile list of sample names
		zgrep -P '^#CHROM\t' ${prefix}Calls_normalized.vcf.gz | cut -f10- | tr "\t" "\n" > ${prefix}Samples.txt
	#generate 1 consensus for each fishSample
		cat ${prefix}Samples.txt | parallel -j $THREADS --no-notice "bcftools consensus -s {} -M '@' -f $ref $prefix''Calls_normalized.vcf.gz -o {}_Consensus.fasta"
		#after running previous line, it is very important to check for consensus seqs that did not assemble properly, they will be short
}


#function to make consensus mtGenome for each group, out and norm
getConsensusByGroup(){
	#assign arguments to variables
		local ref=reference${1}.fasta
		local outCatFILE=outCat${1}-FLTR_F4.bam
		local normCatFILE=normCat${1}-FLTR_F4.bam
		local prefix=$3
		
	#make vcf file
		bcftools mpileup --threads $THREADS -d 20000 -q 30 -Q 20 -A -O z -o out_${prefix}Pile.vcf.gz -f $ref $outCatFILE
		bcftools mpileup --threads $THREADS -d 20000 -q 30 -Q 20 -A -O z -o norm_${prefix}Pile.vcf.gz -f $ref $normCatFILE
	#call genotypes in vcf, force ploidy=haploid
		bcftools call --threads $THREADS -m --ploidy 1 -O z -o out_${prefix}Calls.vcf.gz out_${prefix}Pile.vcf.gz
		bcftools call --threads $THREADS -m --ploidy 1 -O z -o norm_${prefix}Calls.vcf.gz norm_${prefix}Pile.vcf.gz
	#combine snps and indels into 1 multiallelic call and index
		bcftools norm -f $ref -m +any -O z -o out_${prefix}Calls_normalized.vcf.gz out_${prefix}Calls.vcf.gz
		tabix out_${prefix}Calls_normalized.vcf.gz
		bcftools norm -f $ref -m +any -O z -o norm_${prefix}Calls_normalized.vcf.gz norm_${prefix}Calls.vcf.gz
		tabix norm_${prefix}Calls_normalized.vcf.gz
	#compile list of sample names
		zgrep -P '^#CHROM\t' out_${prefix}Calls_normalized.vcf.gz | cut -f10- | tr "\t" "\n" > out_${prefix}Samples.txt
		zgrep -P '^#CHROM\t' norm_${prefix}Calls_normalized.vcf.gz | cut -f10- | tr "\t" "\n" > norm_${prefix}Samples.txt
	#generate 1 consensus for each group
		bcftools consensus -M '@' -f $ref out_$prefix''Calls_normalized.vcf.gz -o out_Consensus.fasta &
		bcftools consensus -M '@' -f $ref norm_$prefix''Calls_normalized.vcf.gz -o norm_Consensus.fasta &
		wait
		#after running previous line, it is very important to check for consensus seqs that did not assemble properly, they will be short
}


#function to extract and align locus to compare between RAD, mtGenomes, and NCBI nucleotide sequences
alignBySample(){
	#assign arguments to variables
		local ef=reference${1}.fasta
		local catfile2=cat${1}-FLTR_F$2.bam
		local prefix=$3
		local positions=$4
		local locus=$5
		eval mtgenomes="$6"

	#get locus from all individuals and align, $positions determines the locus
		cat ${prefix}Samples.txt | parallel -j $THREADS -k --no-notice "echo \>{} && tail -n +2 {}_Consensus.fasta | tr '\n' '\t' | sed 's/\t//g' | sed 's/ */\t/g' | cut -f $positions | sed 's/\t//g' " > ${prefix}RAD_$locus.fasta
	#remove individuals with no sequences or all N
		sed 's/^@@*$//g' ${prefix}RAD_$locus.fasta | grep -P -B 1 '^[ACTGN@]+$' | grep -v '\-\-' > ${prefix}RAD_${locus}_clean.fasta

	#Make fasta from mtGenome sequences
		ls $mtGENOMEs | sed 's/_mtGenome.fasta//g' | parallel -j $THREADS -k --no-notice "echo \>{} && tail -n +2 {}_mtGenome.fasta | tr '\n' '\t' | sed 's/\t//g' | sed 's/ */\t/g' | cut -f $positions | sed 's/\t//g' " > ${prefix}MtGenomes_$locus.fasta

	#Make fasta from other NCBI nucleotide sequences
		#need to code, right now, just download as 1 big fasta
		
	#Combine fastas, make sure to add genbanks
		cat ${prefix}MtGenomes_$locus.fasta ${prefix}RAD_${locus}_clean.fasta > ${prefix}ALL_$locus.fasta 
	#Align all_COI.fasta
		clustalw -infile=${prefix}ALL_$locus.fasta -align -type=DNA -output=NEXUS -outfile=${prefix}ALL_$locus.nex
}

alignByGroup(){
	#assign arguments to variables
		local ef=reference${1}.fasta
		local catfile2=cat${1}-FLTR_F$2.bam
		local prefix=$3
		local positions=$4
		local locus=$5
		eval mtgenomes="$6"

	#get locus from all individuals and align, $positions determines the locus
		echo -e norm\\nout | parallel -j $THREADS -k --no-notice "echo \>{}_RAD_Consensus && tail -n +2 {}_Consensus.fasta | tr '\n' '\t' | sed 's/\t//g' | sed 's/ */\t/g' | cut -f $positions | sed 's/\t//g' " > grp_${prefix}RAD_$locus.fasta
	#remove individuals with no sequences or all N
		sed 's/^@@*$//g' grp_${prefix}RAD_$locus.fasta | grep -P -B 1 '^[ACTGN@]+$' | grep -v '\-\-' > grp_${prefix}RAD_${locus}_clean.fasta

	#Make fasta from mtGenome sequences
		ls $mtGENOMEs | sed 's/_mtGenome.fasta//g' | parallel -j $THREADS -k --no-notice "echo \>{} && tail -n +2 {}_mtGenome.fasta | tr '\n' '\t' | sed 's/\t//g' | sed 's/ */\t/g' | cut -f $positions | sed 's/\t//g' " > ${prefix}MtGenomes_$locus.fasta

	#Make fasta from other NCBI nucleotide sequences
		#need to code, right now, just download as 1 big fasta
		
	#Combine fastas, make sure to add genbanks
		cat ${prefix}MtGenomes_$locus.fasta genbank_COI.fasta grp_${prefix}RAD_$locus.fasta > grp_${prefix}ALL_$locus.fasta 
	#Align all_.fasta
		clustalw -infile=grp_${prefix}ALL_$locus.fasta -align -type=DNA -output=NEXUS -outfile=grp_${prefix}ALL_$locus.nex
}




