rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(seqinr)
library(stringr)

#read in arguments and save into variables
# commandArgs(TRUE)
# args <- commandArgs(TRUE)
inFile <- "fishALL_masked_aligned_clean_eye_mtDNA2_1.fasta"
prefix <- "fishALL_masked_aligned_clean_eye_mtDNA2_1"
#outFile <- args[2]
# num4call <- as.integer(args[3])
# consName <- args[4]
pctMissCall <- 50

#read in data
alignment <- read.alignment(file = inFile, format = "fasta")

#calculate variable values
bp <- nchar(alignment$seq[1])
numMissCall <- round(bp * pctMissCall/100)
missingCalls <- str_count(alignment$seq,pattern="[-Nn]")
pctMissingCalls <- 100*str_count(alignment$seq,pattern="[-Nn]")/bp

#output histograms of missing data (nucleotide calls)
#pdf(file="missingCalls.pdf")
  hist(missingCalls)
  hist(pctMissingCalls)

#make list of seq names and seqs that pass missing data filter
keeperSeqNames <- as.list(alignment$nam[missingCalls < numMissCall])
keeperSeqs <- as.list(alignment$seq[missingCalls < numMissCall])

#convert seqs passing filter into matrix with 1 row for each nucleotide
nuc <- do.call(cbind, strsplit(as.character(keeperSeqs), split=character(0)))

#calculate number of missing calls at each site and make histogram
numN <- sapply(1:dim(nuc)[1], function(x) sum(nuc[x,] == "-" | nuc[x,] == "[nN]"))
hist(numN, xlab="Number of Missing Calls Per Site")

#filter all sites with missing nucleotide calls, then convert from matrix back to list of sequences
keeperNucs <- nuc[numN == 0,]
keeperSeqNucs <- as.list(apply(keeperNucs, 2, function(x) paste(x,collapse="")))

#save filtered alignment as fasta
write.fasta(keeperSeqNucs,keeperSeqNames, file.out = paste(prefix, pctMissCall, ".fasta", sep="_"))

#dev.off()
