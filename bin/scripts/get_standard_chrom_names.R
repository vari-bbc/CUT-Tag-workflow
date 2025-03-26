suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(Rsamtools))

# params from workflow
ref_fasta <- snakemake@params[["ref_fasta"]]
outfile <- snakemake@output[[1]]

# make GenomicRanges from the genome
ref_gr <- as(seqinfo(FaFile(ref_fasta)), "GRanges")

# output the standard chromosome names as a space-separated text file
standardChromosomes(ref_gr) %>% paste(., collapse=" ") %>% writeLines(., outfile)
