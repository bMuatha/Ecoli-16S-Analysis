install.packages("dplyr")
library(dplyr)
#installing required R packages
install.packages("rentrez")
install.packages("seqinr")
install.packages("tidyverse")

BiocManager::install("Biostrings")
BiocManager::install("msa")
aBiocManager::install("ape")
library(rentrez)

#Search and Download Sequences from NCBI
# Search for E. coli 16S rRNA sequences
search <- entrez_search(db="nucleotide",
                        term="Escherichia coli[Organism] AND 16S[Title]",
                        retmax=5)

search$ids   # View IDs

# Download the sequences
seqs <- entrez_fetch(db="nucleotide",
                     id=search$ids,
                     rettype="fasta")

# Save locally
write(seqs, "ecoli_16S_sequences.fasta")

#Read FASTA into R
library(Biostrings)
dna <- readDNAStringSet("ecoli_16S_sequences.fasta")
dna

#Multiple Sequence Alignment (MSA)
library(msa)

alignment <- msa(dna, method="ClustalW")
alignment
msaPrettyPrint(alignment,
               output="pdf",
               showNames="left",
               showLogo="top",
               askForOverwrite=FALSE,
               file="alignment_output.pdf")
png("alignment_output.png", width=1600, height=1200)
msa::msaPrettyPrint(alignment, output="asis")

#GC Content Calculation
library(seqinr)

gc_values <- sapply(as.list(dna), function(x) GC(s2c(as.character(x))))
gc_values

#Build a Phylogenetic Tree
library(ape)

aligned <- msaConvert(alignment, type="ape::DNAbin")
tree <- nj(dist.dna(aligned))
plot(tree, main="Phylogenetic Tree of E. coli 16S Sequences")
