### Detección de SNPs en gen HBB ###
## hemoglobin pseudogen GeneID: 3044
## Length 1607 (from 5241954 to 5243592)

# Librerias:
library(Biostrings)
library(seqinr)

seq_ref <- readDNAStringSet("data/rs10128556.fasta")
seq_ref
class(seq_ref)

# Crear función para converti el objeto de clase dnastringset en un dataframe:
dss2df <- function(dss){
  data.frame(width=width(dss),
             seq=as.character(dss), 
             names=names(dss))
}

# Convertir el objeto:
seq2 <- dss2df(seq_ref) 
print(seq2$seq)

seq3 <- seq2$seq
dim(seq3)

seq_letters <- unlist(strsplit(seq3,""))
seq_letters
length(seq_letters)
class(seq_letters)

# Crear archivo fasta de la secuecnia que contiene el snp:

snp_txt <- c("AGTTTGAGCAACTCTCACCATTATGGGCTAGATATCAGCTTCAATATGACTGGCAGCCCTGCACCTCCCA
TTGTGTCCTATCCTTCTCTTACTTGCTATGTCAACTCACTACCCCAACATATTGTGATTGTTTCATTTTT
TTTTTAGAGACTTCGTCCTTTTAAAAAGTAGGATAATAGTACTTCAAGGAACAGTAATGGA")
write.fasta(sequences = snp_txt, names = names("rs10128556"), 
            nbchar = 60, file.out = "rs10128556.fasta")

# Alineamiento de la secuencia que contiene el snp con el gen:
seq_snp <- readDNAStringSet("data/rs10128556.fasta")

sigma <- nucleotideSubstitutionMatrix(match = 2, 
                                      mismatch = -1, 
                                      baseOnly = TRUE)

alignment <- pairwiseAlignment(seq_ref, seq_snp, 
                               substitutionMatrix = sigma, 
                               gapOpening -2, 
                               gapExtension = -8, 
                               scoreOnly = FALSE)
print(alignment)

seq1.1 <- readDNAStringSet("data/1.1.fasta")
alignment_seq1.1 <- pairwiseAlignment(seq1.1, seq_snp, 
                                      substitutionMatrix = sigma, 
                                      gapOpening -2, 
                                      gapExtension = -8, 
                                      scoreOnly = FALSE)
alignment_seq1.1

seq2.1 <- readDNAStringSet("data/2.1.fasta")
alignment_seq2.1 <- pairwiseAlignment(seq2.1, seq_snp, 
                                      substitutionMatrix = sigma, 
                                      gapOpening -2, 
                                      gapExtension = -8, 
                                      scoreOnly = FALSE)
alignment_seq2.1
