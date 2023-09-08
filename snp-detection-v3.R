# Detección de SNPs descritos en literatura a partir de los datos tratados
# Se ingresará una matrix que corresponde al resultado del alineamientom y corte
# de la secuencia de interés sobre la socuencia de referencia de NCBI (NC_000011)
# se probará usando el SNP rs10128556 (C>T)

rm(list=ls())

# Load packages:
library(Biostrings)
library(DECIPHER)
library(adegenet)
library(ips)
library(seqinr)

# Read matrix
seqs <- read.table("data/processed/test1.txt", header = TRUE, row.names = 1)
head(seqs)
sub_as_string <- seqs[,1]
seq_string <- paste(sub_as_string, collapse = " ")

snp <- read.fasta(file = "data/raw/rs10128556_snp.fa", as.string=TRUE)
snq_seq <- c("tatccttctcttacttgctatgtcaactcactaccccaacatattgtg")

result <- sapply(seq_string, patter = snq_seq, overlap=TRUE)

sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1,
                                      baseOnly = TRUE)

result <- pairwiseAlignment(pattern = seq_sstring,
                            subject = snp,
                            substitutionMatrix = sigma, 
                            gapOpening -2, 
                            gapExtension = -8, 
                            scoreOnly = FALSE)
