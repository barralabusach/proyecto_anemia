## En este script se alinea el gen de BCL11A con uno de los SNPs descritos.

library(Biostrings)
library(ape)
library(pegas)
library(msa)
library(adegenet)
library(EnvNJ)
library(seqinr)

bcl11a <- readDNAStringSet("data/bcl11a.fasta")
snp1 <- readDNAStringSet("data/rs7581162 .fa")
bcl11a <- read.fasta("data/bcl11a.fasta")
snp1 <- read.fasta("data/rs7581162.fa")

snp1 <- snp1[[1]]
write.fasta(name="rs7581162", sequences=snp1, file.out="rs7581162.fasta")


sequences <- fastaconc(otus = c("bcl11a", "rs7581162"),
                       inputdir = "data/.",
                       out.file = "bcl11a-snp1.fasta")

sequences <- readDNAStringSet("data/bcl11a-snp1.fasta")
bcl11a.alignment <- msa(sequences)
print(bcl11a.alignment, show="complete")

bcl11a <- bcl11a[[1]]
length(bcl11a)
snp1 <- snp1[[1]]
length(snp1)
seq_matriz <- matrix(bcl11a)
dim(seq_matriz)

library(stringr)
