rm(list=ls())

# librerias
library(bioseq)
library(tidyverse)
library(EnvNJ)

# crear fasta con todas las secuencias
sequences <- EnvNJ::fastaconc(
  otus = c("sample1", "sample2", "sample3",
           "sample4", "sample5", "sample6", 
           "sample7", "sample8", "sample9", 
           "sample10", "sample11"),
  inputdir = "data/raw/", 
  out.file  = "data/interm/conc_samples.fasta"
  )

# leer archivo fasta
seq_data <- read_fasta("data/interm/conc_samples.fasta")
seq_data

seq_nchar(seq_data) %>%
  range()

samples_tab <- tibble(sequence = seq_data)

samples_tab$names <- names(c("sample1", "sample2", "sample3",
                            "sample4", "sample5", "sample6", 
                            "sample7", "sample8", "sample9", 
                            "sample10", "sample11"))


