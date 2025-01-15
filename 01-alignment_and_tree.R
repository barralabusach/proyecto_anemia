## Script para realizar clasificación de las secuencias de hemoglobina de una muestra de población chilena ##

# Clear memory space:
rm(list = ls())

### Install Bioconductor to install the needed packages:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

# load packages:
if(!require("Biostrings")){BiocManager::install("Biostrings")}
if(!require("msa")){BiocManager::install("msa")}
if(!require("EnvNJ")){BiocManager::install("EnvNJ")}
if(!require("ggtree")){BiocManager::install("ggtree", force = TRUE)}
if(!require("ggplot2")){BiocManager::install("ggplot2")}
if(!require("adegenet")){BiocManager::install("adegenet")}
if(!require("ape")){BiocManager::install("ape")}
if(!require("ips")){BiocManager::install("ips")}
if(!require("bios2mds")){BiocManager::install("bios2mds")}

# Entries:
# Create one file with several sequences:
sequences <- EnvNJ::fastaconc(
  otus = c(
    "sample1-hbb", "sample2-hbb", "sample3-hbb", "sample4-hbb", "sample5-hbb",
    "sample6-hbb", "sample7-hbb", "sample8-hbb", "sample9-hbb", "sample10-hbb",
    "sample11-hbb"
  ),
  inputdir = "data/",
  out.file  = "data/concatenates_sequences.fasta"
)

# Read file:
sequences <- Biostrings::readDNAStringSet(
  "data/interm/concatenates_sequences.fasta"
)

# Multiple sequence alignment CluscalW por default:
seq_alignment <- msa(sequences, method = "ClustalW")
print(seq_alignment, show = "complete")

# Convert alignment file to fasta file and export:
fasta_file <- msaConvert(seq_alignment, type = c("bios2mds::align"))
export.fasta(
  fasta_file,
  outfile = "data/processed/outfile.fas",
  ncol(seq_alignment),
  open = "w"
)

# Read aligned data from converted file above, convert to DNAbin object:
aligned_data <- adegenet::fasta2DNAbin("data/processed/outfile.fas")

# After align, trim the sequence:
trimed_seq <- ips::trimEnds(aligned_data) # de nbin a trimed_seq

# Convert DNAbin to alignment format:
align_format <- ape::as.alignment(trimed_seq) # de an a align_format

# Convert alignment to matrix:
align_matrix <- as.matrix(align_format) # de nm a align_matrix

# Extraction of the sample names:
sample_names <- as.matrix(labels(aligned_data))
sample_names

# Computing distance by ape package with K80 model derived by Kimura (1980):
# Kimura model assumes iqual base frequencies and accounts for the difference
# between transitions and transversions with one parameter.
comp_dist <- ape::dist.dna(trimed_seq, model = "K80")
tree <- nj(comp_dist)

ggt <- ggtree(
  tree, 
  cex = 0.8,
  aes(color = branch.length)
  ) +
  scale_color_continuous(
    high = "lightskyblue1", 
    low = "red"
    ) +
  geom_tiplab(
    align = TRUE,
    size = 2
    ) + 
  geom_treescale(
    y = 0.5,
    color = "coral4",
    fontsize = 4
    )

ggt <- ggtree(
  tree, 
  branch.length = "none",
  aes(color = branch.length)
) +
  scale_color_continuous(
    high = "lightskyblue1", 
    low = "red"
  ) +
  geom_tiplab(
    align = TRUE,
    size = 2
  ) + 
  geom_treescale(
    y = - 5,
    color = "coral4",
    fontsize = 4
  )

tree_plot <- msaplot(
  ggt,
  trimed_seq,
  offset = 2,
  width = 1,
  height = 0.5,
  color = c(
    rep("rosybrown", 1),
    rep("sienna1", 1),
    rep("lightgoldenrod1", 1),
    rep("lightskyblue1", 1),
    rep("green", 1),
    rep("yellow", 1)
  )
)
print(tree_plot)
