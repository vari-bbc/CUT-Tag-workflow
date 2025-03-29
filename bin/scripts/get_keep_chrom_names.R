suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(stringr))

# params from workflow
ref_fasta <- snakemake@input[["ref_fasta"]]
outfile <- snakemake@output[[1]]
keep_std_chroms <- snakemake@params[['keep_std_chroms']]  #logical
rm_chroms <- snakemake@params[["rm_chroms"]] # comma-separated string or "foobar" if no chroms indicated for removal

# make GenomicRanges from the genome
ref_gr <- as(seqinfo(FaFile(ref_fasta)), "GRanges")
if(keep_std_chroms){
    ref_gr <- keepStandardChromosomes(ref_gr, pruning.mode="coarse")
}
if(rm_chroms != "foobar"){
    rm_chroms <- str_split_1(rm_chroms, ",")
    ref_gr <- dropSeqlevels(ref_gr, rm_chroms, pruning.mode="coarse")
}


# output the standard chromosome names as a space-separated text file
seqlevels(ref_gr) %>% paste(., collapse=" ") %>% writeLines(., outfile)
