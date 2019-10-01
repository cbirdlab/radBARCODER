#!/usr/bin/env Rscript
#script to create consensus sequence

#read in arguments and save into variables
	commandArgs(TRUE)
	args <- commandArgs(TRUE)
		inFile <- args[1]
		outFile <- args[2]
		num4call <- as.integer(args[3])
		consName <- args[4]

library(seqinr)

#read in data
	alignment <- read.alignment(file = inFile, format = "fasta")

#make consensus profile, the profile option will allow me to control how the consensus is made
	alignment_profile <- consensus(alignment, method = "profile", warn.non.IUPAC = FALSE, type = "DNA")
	numBP <- length(alignment_profile[1,])

#make function to call consensus nucleotide
	call_consensus_nucleotide <- function(nuc){
		#load variables
		ACTG <- nuc[names(nuc) != "-" & names(nuc) != "n"]
		numACTG <- sum(ACTG)

		#if there are at least num4call good nucleotide calls, then call the nucleotide, else call n
		if(numACTG >= num4call){
			nuc_call <- toupper(names(ACTG[which.max(ACTG)]))
		} else {
			nuc_call <- "-"
		}
		return(nuc_call)
	}

#if there are at least 2 calls that are not N then call the base
	alignment_consensus <- sapply(1:numBP, FUN=function(x) call_consensus_nucleotide(alignment_profile[,x]))

#save consensus to file
	write.fasta(alignment_consensus, names = consName, file.out = outFile, open = "w")