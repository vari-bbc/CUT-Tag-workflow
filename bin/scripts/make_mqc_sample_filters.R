##########################################
### Variables set by snakemake pipeline###
##########################################

samplesheet <- snakemake@input[["samplesheet"]]
out_sample_filts <- snakemake@output[["sample_filts"]]
noncontrol_samples <- snakemake@params[["noncontrol_samples"]]
##########
##########
##########

# load packages
library(dplyr)
library(stringr)
library(readr)

meta <- read_tsv(samplesheet) %>%
    dplyr::mutate(mqc_group = case_when(
        sample %in% str_split_1(noncontrol_samples, ",") ~ "Controls",
        .default = sample_group
))
    dplyr::select(-fq1, -fq2) %>%
    unique()

fileConn <- file(out_sample_filts)
mqc_groups <- unique(meta$mqc_group)
for (i in 1:length(mqc_groups)){
    curr_group <- mqc_groups[i]
    curr_samples <- meta$sample[which(meta$mqc_group == curr_group), ]
    outline <- paste(c(curr_group, "show", curr_samples), sep="\\t")
    writeLines(outline, fileConn)
}

close(fileConn)

sessionInfo()
