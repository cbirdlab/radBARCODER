#!/usr/bin/env Rscript

#rm(list=ls())

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(seqinr)
library(stringr)

#read in arguments and save into variables
	commandArgs(TRUE)
		args <- commandArgs(TRUE)
		inFile <- paste(args[1],".fasta",sep="")
		outFile <- paste(args[1],"_",args[2],".fasta", sep="")
		pctMissCall <- as.numeric(args[2])

#read in data
	alignment <- read.alignment(file = inFile, format = "fasta")

#calculate variable values
	bp <- nchar(alignment$seq[1])
	numMissCall_bp <- round(bp * pctMissCall/100)
	missingCalls_bp <- str_count(alignment$seq,pattern="[-Nn]")
	pctmissingCalls_bp <- 100*str_count(alignment$seq,pattern="[-Nn]")/bp

#output histograms of missing data (nucleotide calls)
pdf(file=paste("missingCalls_", pctMissCall, ".pdf", sep=""))
  hist(missingCalls_bp, 
       main="Number Ambiguous/Missing/Indel Nucleotide Calls Per Genome",
       sub = paste("Alignment =", bp, "bp", sep=" "),
       ylab="Number of Genomes", 
       xlab="# Missing/Indel/Ambiguous Nucleotide Calls")
  hist(pctmissingCalls_bp, 
       main=paste(length(keeperSeqs), "Genomes With < ", pctMissCall, "% Missing/Indel/Ambiguous Nucleotide Calls", sep=""),
       ylab="Number of Genomes", 
       xlab="% Missing/Indel/Ambiguous Nucleotide Calls")

#make list of seq names and seqs that pass missing data filter
	keeperSeqNames <- as.list(alignment$nam[missingCalls_bp < numMissCall_bp])
	keeperSeqs <- as.list(alignment$seq[missingCalls_bp < numMissCall_bp])

#convert seqs passing filter into matrix with 1 row for each nucleotide
	nuc <- do.call(cbind, strsplit(as.character(keeperSeqs), split=character(0)))

#calculate number of missing calls at each site and make histogram
	numN <- sapply(1:dim(nuc)[1], function(x) sum(nuc[x,] == "-" | nuc[x,] == "[nN]"))
	keeperSites <- which(numN == 0)
	keeperSites <- split(keeperSites, cumsum(c(1, diff(keeperSites) != 1)))
	keeperSiteRanges <- lapply(keeperSites, range)
	write.csv(keeperSiteRanges, file=paste("genomePositionsRemaining_",pctMissCall,".csv", sep=""), row.names=FALSE)
	bp_remaining <- as.integer(table(numN)[1])
	
	plot(numN,
	     main=paste(bp_remaining, "Nucleotides With No Missing/Ambiguous/Indel Calls in", length(keeperSeqs), "Genomes", sep=" "),
	     xlab="Position in Genome",
	     ylab="Number of Missing/Indel/Ambiguous Calls Per Site")
	hist(numN, 
	     main="Number of Missing/Indel/Ambiguous Calls Per Site",
	     sub=paste(length(keeperSeqs), "Genomes Remain", sep=" "),
	     xlab="# of Missing/Indel/Ambiguous Calls Per Site", 
	     ylab="Number of Sites", 
	     breaks=max(numN))

#filter all sites with missing nucleotide calls, then convert from matrix back to list of sequences
	keeperNucs <- nuc[numN == 0,]
	keeperSeqNucs <- as.list(apply(keeperNucs, 2, function(x) paste(x,collapse="")))

#save filtered alignment as fasta
	write.fasta(keeperSeqNucs,keeperSeqNames, file.out = outFile)

dev.off()
