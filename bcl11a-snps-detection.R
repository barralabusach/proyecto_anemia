# Script para la detección de SNPs en la secuencia del gen de BCL11A
# Se usará como molde la secuencia concenso

rm(list=ls())

# Librerias:
library(bioseq)
library(tidyverse)
library(EnvNJ)

# Leer archivo de secuencia de gen bcl11a:
bcl11a <- bioseq::read_fasta("data/raw/bcl11a-template.fa")
bcl11a

bcl11a_tab <- dplyr::tibble(label = names(bcl11a), sequence = bcl11a)

x_dna <- dna(bcl11a_tab$sequence)
x_dna
is_dna(x_dna)

# Flancos de SNP rs4671393 de 25nt:
# CCAGTGCTGTGGACAGCAAAGCTTCGTGCAGGAAATTAAGATTCCCCCTG
templ <- dna("CCAGTGCTGTGGACAGCAAAGCTTCGTGCAGGAAATTAAGATTCCCCCTG")
x <- seq_complement(templ)

seq_detect_pattern(x_dna, dna("CCAGTGCTGTGGACAGCAAAGCTTCGTGCAGGAAATTAAGATTCCCCCTG"))
seq_detect_pattern(x_dna, x)

bcl11a_tab <- bcl11a_tab %>%
  mutate(SNP = seq_detect_pattern(x_dna, x))



