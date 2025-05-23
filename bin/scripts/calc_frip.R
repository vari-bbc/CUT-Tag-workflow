##########################################
### Variables set by snakemake pipeline###
##########################################

samplesheet <- snakemake@input[["samplesheet"]]
frag_count_files <- snakemake@input[["frag_counts"]]
out_all_frip <- snakemake@output[["all_merge"]]
out_group_frip <- snakemake@output[["group_merge"]]
peak_type <- snakemake@params[['peak_type']]
noncontrol_samples <- snakemake@params[["noncontrol_samples"]]
##########
##########
##########

# load packages
library(dplyr)
library(stringr)
library(readr)

meta <- read_tsv(samplesheet) %>%
    dplyr::filter(sample %in% str_split_1(noncontrol_samples, ",")) %>%
    dplyr::arrange(enriched_factor, sample_group, sample) %>%
    dplyr::select(sample, sample_group) %>%
    unique()

stopifnot(length(meta$sample) == length(unique(meta$sample)))

by_group_frip <- lapply(setNames(nm=meta$sample), function(x){

    infile <- file.path("analysis/frags_over_peaks/", peak_type, meta$sample_group[which(meta$sample == x)], paste0(x, ".txt"))
    print(infile)
    stopifnot(str_replace_all(infile, "\\/+", "/") %in% str_replace_all(frag_count_files, "\\/+", "/"))

    message(str_glue("By-group FRiP for {x} from {infile}"))
    
    lines <- readLines(infile)
    lines <- unlist(lapply(lines, as.numeric))
    
    data.frame(`Sample Name` = x, peak_frags = lines[2], nonpeak_frags = lines[1] - lines[2])
})

do.call(rbind, by_group_frip) %>% 
    write_csv(., out_group_frip)


all_frip <- lapply(setNames(nm=meta$sample), function(x){
    infile <- file.path("analysis/frags_over_peaks/", peak_type, "all", paste0(x, ".txt"))

    stopifnot(str_replace_all(infile, "\\/+", "/") %in% str_replace_all(frag_count_files, "\\/+", "/"))

    message(str_glue("All-samples FRiP for {x} from {infile}"))
    
    lines <- readLines(infile)
    lines <- unlist(lapply(lines, as.numeric))

    data.frame(`Sample Name` = x, peak_frags = lines[2], nonpeak_frags = lines[1] - lines[2])
})

do.call(rbind, all_frip) %>% 
    write_csv(., out_all_frip)

sessionInfo()
