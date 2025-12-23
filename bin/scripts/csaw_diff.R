options(warn=1) # show all warnings

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

contrasts_file <- snakemake@input[['contrasts']]
filt <- snakemake@input[['filt_small_wins']]
binned <- snakemake@input[['bin_counts']]
filt.spikein <- snakemake@input[['filt_small_wins_spikein']]
binned.spikein <- snakemake@input[['bin_counts_spikein']]

peaks <- snakemake@input[['merged_peaks']]
norm_type <- snakemake@wildcards[['norm_type']]
enriched_factor <- snakemake@wildcards[['enriched_factor']]

# output RDS objects
results_rds <- snakemake@output[['results_rds']]
results_SEs_rds <- snakemake@output[['results_se_rds']]
edgeR_rds <- snakemake@output[['edgeR_rds']]
peaks_res_rds <- snakemake@output[['peaks_res_rds']]
combined_rds <- snakemake@output[['combined_rds']]
peaks_anno_rds <- snakemake@output[['peaks_anno_rds']]

orgdb_pkg <- snakemake@params[['orgdb']]
txdb_pkg <- snakemake@params[['txdb']]
promo_start <- as.integer(snakemake@params[['promo_start']])
promo_end <- as.integer(snakemake@params[['promo_end']])


##########
##########
##########

# load packages
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(csaw))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(ChIPseeker))
if (!require(orgdb_pkg, quietly = TRUE))
    BiocManager::install(orgdb_pkg)
if (!require(txdb_pkg, quietly = TRUE))
    BiocManager::install(txdb_pkg)
suppressPackageStartupMessages(library(orgdb_pkg, character.only = TRUE))
suppressPackageStartupMessages(library(txdb_pkg, character.only = TRUE))

filt <- readRDS(filt)
#binned <- readRDS(binned)

contrasts <- read_tsv(contrasts_file)
contrasts <- contrasts %>%
    dplyr::filter(group1 %in% filt$sample_group & group2 %in% filt$sample_group) %>%
    unique()

if (nrow(contrasts) < 1){
    stop(str_glue("No valid contrasts specified for {enriched_factor}"))
}

peaks <- import(peaks)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
message("Range of peak widths is ", paste(range(width(peaks)), collapse="-"))

# annotate all peaks with ChIPSeeker
chipseeker_args <- list(peak=peaks, annoDb = orgdb_pkg, TxDb = eval(as.name(txdb_pkg)), tssRegion = c(promo_start, promo_end))
chipseeker_args
peaks <- as.GRanges(do.call(annotatePeak, chipseeker_args))
write_rds(peaks, peaks_anno_rds)

if (norm_type == "bkgd"){

    binned <- readRDS(binned)
    stopifnot(identical(colnames(filt), colnames(binned)))
    filt$norm.factors <- binned$norm.factors

} else if (norm_type == "bkgd.spikein"){

    binned.spikein <- readRDS(binned.spikein)
    stopifnot(identical(colnames(filt), colnames(binned.spikein)))
    filt$norm.factors <- binned.spikein$norm.factors

} else if (norm_type == "hiAbund.spikein"){

    filt.spikein <- readRDS(filt.spikein)
    stopifnot(identical(colnames(filt), colnames(filt.spikein)))
    filt$norm.factors <- filt.spikein$norm.factors

} else if (norm_type == "hiAbund"){

    # do nothing; the norm_factors in 'filt' are already high abundance window normalization factors
    stopifnot("norm.factors" %in% colnames(colData(filt)))
    
} else {
    stop("Invalid norm_type")
}

y <- asDGEList(filt, samples=colData(filt)[, c("sample","sample_group","enriched_factor")])
colnames(y) <- y$samples$sample

design <- model.matrix(~0+sample_group, data = y$samples)
colnames(design) <- str_remove(colnames(design), "^sample_group")

# vector of strings describing contrast in the from "B-A" where B and A are levels in the sample_group column
comps <- paste(contrasts[["group2"]], contrasts[["group1"]], sep="-")
#print(design)
#print(contrasts)
stopifnot(all(c(contrasts[["group1"]], contrasts[["group2"]]) %in% colnames(design)))
cons <- makeContrasts(contrasts=comps, levels=design)

y <- estimateDisp(y, design)
summary(y$trended.dispersion)

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.post)

pdf(file.path(figs_dir, "squeezing.pdf"), width=9, height=5)
par(mfrow=c(1,2))
o <- order(y$AveLogCPM)
plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
     ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
     ylab=("Biological coefficient of variation"))
plotQLDisp(fit)
dev.off()

results <- lapply(setNames(nm=colnames(cons)), function(x){
  glmQLFTest(fit, contrast = cons[, x])
})

results_SEs <- lapply(results, function(x) {
  se <- filt
  rowData(se) <- cbind(rowData(se), x$table)
  se
})

saveRDS(results, results_rds)
saveRDS(results_SEs, results_SEs_rds)
saveRDS(list(y=y, fit=fit), edgeR_rds)

peaks_res <- lapply(results_SEs, function(x) {
  out <- overlapResults(x, regions=peaks, tab=rowData(x))
  out$best$rep.logCPM.manual <- rowData(x)$logCPM[out$best$rep.test]
  out$combined$rep.logCPM.manual <- rowData(x)$logCPM[out$combined$rep.test]
  
  out
})

combined <- lapply(peaks_res, function(x) {
  gr <- x$regions
  mcols(gr) <- cbind(mcols(gr), x$combined)
  gr
})

saveRDS(peaks_res, peaks_res_rds)
saveRDS(combined, combined_rds)

# session info
sessionInfo()
