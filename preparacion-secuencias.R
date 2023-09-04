## Script para tratar las secuencias usando com referencia la secuencia NC_000011.10

# Set env:
rm(list=ls())
setwd("C:/Users/cesar/OneDrive/Proyectos/Anemia")

# Load packages:
library(Biostrings)
library(DECIPHER)
library(adegenet)
library(ips)

pattern_seq <- readDNAStringSet("data/raw/NC_000011.10.fasta") # pattern
subject_seq <- readDNAStringSet("data/raw/sample1.fasta") # subject

# Nedleman-Wunsch alignment alogrithm:
sigma <- nucleotideSubstitutionMatrix(
      match = 2,
      mismatch = -1,
      baseOnly = TRUE
      )

seq_alignment <- pairwiseAlignment(
      pattern = pattern_seq, 
      subject =  subject_seq,
      substitutionMatrix = sigma, 
      gapOpening -2, 
      gapExtension = -8, 
      scoreOnly = FALSE
      )

# Extract the alignet pattern and show in browser:
sequences_char <- c(
      aligned(pattern(seq_alignment)), 
      aligned(subject(seq_alignment))
      )

#as.character(sequences_char)
#class(sequences_char)
BrowseSeqs(sequences_char)

# Convert alignment in fasta file:
. <- c(
      as(subject(seq_alignment), "DNAStringSet"), 
      as(pattern(seq_alignment), "DNAStringSet")
      )
names(.) = c("subject", "pattern")

writeXStringSet(., "data/interm/PairwiseAlignment.fasta")

# Trim sequences:
aligned_data <- adegenet::fasta2DNAbin("data/interm/PairwiseAlignment.fasta")

trim_seq <- ips::trimEnds(aligned_data)
as_align <- as.alignment(trim_seq)

as_matrix <- as.matrix(as_align)
as_matrix <- t(as_matrix)
dim(as_matrix)
class(as_matrix)
View(as_matrix)

# Save matrix
write.table(as_matrix, file = "data/processed/test1.txt")
