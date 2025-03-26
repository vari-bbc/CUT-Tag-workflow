##########################################
### Variables set by snakemake pipeline###
##########################################

# for debugging
#save.image("/varidata/research/projects/bbc/research/JONR_20220120_cutnrun_VBCS-694/foo.Rdata")
#print(str(snakemake))

figs_dir <- snakemake@output[['figs_dir']]
dir.create(figs_dir, recursive=TRUE)

# cores for parallel processing
bp_cores <- snakemake@threads-1 

# blacklist BED file
blacklist_file <- snakemake@params[["blacklist"]]

# file with the standard chromosome names
# file should contain only 1 line, which is space delimited
std_chroms_file <- snakemake@input[["std_chroms"]]

# read in mitochondrial chromosome name
mito_chr <- snakemake@params[["mito_chr"]]


# bam files to read in
bam.files <- snakemake@input[["bams"]]

names(bam.files) <- snakemake@params[["samp_names"]] 

# output RDS objects
binned_rds <- snakemake@output[['binned']]
small_wins_rds <- snakemake@output[['small_wins']]
filt_small_wins_rds <- snakemake@output[['filt_small_wins']]
global_filt_rds <- snakemake@output[['global_filt']]
outCompNorm <- snakemake@output[["bkgrd_scale"]]
outEffNorm <- snakemake@output[["hiAbund_scale"]]
outCompNormLink <- snakemake@output[["bkgrd_scale_link"]]
outEffNormLink <- snakemake@output[["hiAbund_scale_link"]]

# csaw params
bkgd_bin_width <- 10000
hi_abund_win_width <- snakemake@params[["window_width"]]
win_filter_fold <- 3 # fold change from background to keep as high abundance windows
minq <- 20
dedup <- FALSE
hi_abund_th <- log2(win_filter_fold)

##########
##########
##########

# load packages
suppressPackageStartupMessages(library(csaw))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))

# read in standard chromosomes file as a vector, and remove mitochondrial chromosome name
chroms_keep_vec <- strsplit(readLines(std_chroms_file)[1], " ")[[1]] %>% 
    grep(pattern=paste0("^", mito_chr, "$"), x=., perl=TRUE, value=TRUE, invert=TRUE)

message("Keeping only these chromosomes: ", paste(chroms_keep_vec, collapse=", "))


# import blacklist as genomicranges
blacklist_gr <- import(blacklist_file)

# Set params for reading in BAM files
param <- readParam(minq=minq, dedup=dedup, pe="both", restrict=chroms_keep_vec, discard=blacklist_gr)

# Eliminate composition biases
# TMM on binned counts

binned <- windowCounts(bam.files, bin=TRUE, width=bkgd_bin_width, param=param, BPPARAM=MulticoreParam(bp_cores)) # note that windows with less than 'filter' number of reads summed across libraries are removed
message("For TMM on background bins, ", length(binned), " bins were used.")

binned <- normFactors(binned, se.out=TRUE)
saveRDS(binned, binned_rds) 

# code adapted from https://www.biostars.org/p/394434/
# 'se' is csaw generated ranged summarized experiment with norm.factors and totals in the coldata
# 'outfile' is the output file containing the colData
calc_and_write_final_norm_facs <- function(se, outfile){
    se$final.factors <- ((se$norm.factors * colData(se)$totals) / 1000000)^-1 # Recall that the norm factors from edgeR/csaw must be multiplied by the library size to also correct for library size; this is different from and not need with the DESeq2 norm factors.
    write.table(x = as.data.frame(colData(se)) %>% tibble::rownames_to_column("sample"),
                file = outfile, sep="\t", quote = FALSE,
                col.names = TRUE, row.names = FALSE)

}
## turn off exponential notation to avoid rounding errors
options(scipen=999)

calc_and_write_final_norm_facs(se=binned, outfile=outCompNorm)
system(str_glue("ln -sr {outCompNorm} {outCompNormLink}"))

# Eliminate efficiency biases
# TMM on high abundance regions

small_wins <- windowCounts(bam.files, width=hi_abund_win_width, param=param, BPPARAM=MulticoreParam(bp_cores)) # note that windows with less than 'filter' number of reads summed across libraries are removed
saveRDS(small_wins, small_wins_rds) # for faster debugging

filter_wins <- filterWindowsGlobal(small_wins, binned)
saveRDS(filter_wins, global_filt_rds)

pdf(file.path(figs_dir, "filter_th.pdf"), width=6, height=5)
hist(filter_wins$filter, main="", breaks=50,
    xlab="Log-fold change from global background")
abline(v=hi_abund_th, col="red")
dev.off()

keep <- filter_wins$filter > hi_abund_th # filter for greater than log2('win_filter_fold') fold difference compared to background

keep_num <- sum(keep)
tot_wins <- length(keep)
message("For TMM on high abundance normalization, ", keep_num, " out of ", tot_wins, " (", round(keep_num/tot_wins, 2), ")", " windows were kept.")
filtered.wins <- small_wins[keep, ]

filtered.wins <- normFactors(filtered.wins, se.out=TRUE) 
saveRDS(filtered.wins, filt_small_wins_rds) 

calc_and_write_final_norm_facs(se=filtered.wins, outfile=outEffNorm)
system(str_glue("ln -sr {outEffNorm} {outEffNormLink}"))
