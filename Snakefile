import pandas as pd
import numpy as np
import os
import re
import itertools
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("7.25.0")

##### Editable variables #####

configfile: "bin/config.yaml"

bt2_index = config['ref']['index']

##### load config and sample sheets #####
samplesheet="bin/samples.tsv"
units = pd.read_table(samplesheet, dtype={"sample" : str, "sample_group" : str })
units['se_or_pe'] = ["SE" if x else "PE" for x in units['fq2'].isnull()]

samples = units[["sample","control","sample_group","enriched_factor","se_or_pe"]].drop_duplicates()
if not samples['sample'].is_unique:
    raise Exception('A sample has more than one combination of control, sample_group, enriched_factor, and/or se_or_pe.')

snakemake_dir = os.getcwd() + "/"

# make a tmp directory for analyses
tmp_dir = os.path.join(snakemake_dir, "tmp")
if not os.path.exists(tmp_dir):
    os.mkdir(tmp_dir)


rule all:
    input:
        #expand("analysis/bwamem/{sample.sample}.bam", sample=samples.itertuples()),
        #expand("analysis/filt_bams/{sample.sample}_filt_alns.bam", sample=samples.itertuples()),
        #expand("analysis/bowtie2/{sample.sample}.sorted.bam", sample=samples.itertuples()),
        #expand("analysis/bed_files/{sample.sample}.bed.gz", sample=samples.itertuples())
        expand("analysis/bigwig_files/{sample.sample}.bw", sample=samples.itertuples())
def get_orig_fastq(wildcards):
    if wildcards.read == "R1":
            fastq = expand("raw_data/{fq}", fq = units[units["sample"] == wildcards.sample]["fq1"].values)
    elif wildcards.read == "R2":
            fastq = expand("raw_data/{fq}", fq = units[units["sample"] == wildcards.sample]["fq2"].values)
    return fastq 

rule rename_fastqs:
    """
    Rename fastqs by biologically meaningful name. Concatenate different runs of same library.
    """
    input:
        get_orig_fastq
    output:
        "analysis/renamed_data/{sample}_{read}.fastq.gz"
    benchmark:
        "benchmarks/rename_fastqs/{sample}_{read}.txt"
    params:
        num_input=lambda wildcards, input: len(input),
        input=lambda wildcards, input: ["'" + x + "'" for x in input]
    threads: 1
    resources:
        mem_gb=8,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
    shell:
        """
        if [ {params.num_input} -gt 1 ]
        then
            cat {params.input} > {output}
        else
            ln -sr {params.input} {output}
        fi

        """

rule fastqc:
    """
    Run fastqc on raw_data/ files.
    """
    input:
        "analysis/renamed_data/{fq_pref}.fastq.gz"
    output:
        html="analysis/fastqc/{fq_pref}_fastqc.html",
        zip="analysis/fastqc/{fq_pref}_fastqc.zip"
    params:
        outdir="analysis/fastqc/"
    benchmark:
        "benchmarks/fastqc/{fq_pref}.txt"
    envmodules:
        config['modules']['fastqc']
    threads: 1
    resources:
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        fastqc --outdir {params.outdir} {input}
        """

rule cutadapt:
    input:
       # expand("analysis/renamed_data/{{sample}}_R{read}.fastq.gz", read=["1","2"])
        R1="analysis/renamed_data/{sample}_R1.fastq.gz",
        R2="analysis/renamed_data/{sample}_R2.fastq.gz"
    output:
        R1_trimmed="analysis/trim_fq/{sample}_R1_trimmed.fq.gz",
        R2_trimmed="analysis/trim_fq/{sample}_R2_trimmed.fq.gz"
    benchmark:
        "benchmarks/trim_fq/{sample}.txt"
    envmodules:
        config['modules']['cutadapt']
    threads: 4
    resources:
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        cutadapt -j 8 -m 20 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o {output.R1_trimmed} -p {output.R2_trimmed} {input.R1} {input.R2}
        """

rule bowtie2:
    input:
        R1_trimmed="analysis/trim_fq/{sample}_R1_trimmed.fq.gz",
        R2_trimmed="analysis/trim_fq/{sample}_R2_trimmed.fq.gz"   
    output:
        #sam="analysis/bowtie2/{sample}.sam",
        sorted_bam="analysis/bowtie2/{sample}.sorted.bam",
        unsorted_bam="analysis/bowtie2/{sample}.bam",
        outbai="analysis/bowtie2/{sample}.sorted.bam.bai",
        idxstat="analysis/bowtie2/{sample}.bam.idxstat"
    benchmark:
        "benchmarks/bowtie2/{sample}.txt"
    params:
        bt2_index=bt2_index,
        samblaster_params=lambda wildcards: "--addMateTags" if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "--ignoreUnmated"
    threads: 16
    envmodules:
        config['modules']['bowtie2'],
        config['modules']['samblaster'],
        config['modules']['samtools'],
    resources:
        mem_gb=180,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        bowtie2 -p 16 --very-sensitive-local --soft-clipped-unmapped-tlen --no-mixed --no-discordant --dovetail --phred33 \
        -I 10 -X 1000 -x {params.bt2_index} -1 {input.R1_trimmed} -2 {input.R2_trimmed} | samtools view -bS \
        -@ {threads} \
        -O "BAM" \
        -o {output.unsorted_bam} 
        
        samtools sort {output.unsorted_bam} -m 6G -O "bam" -o {output.sorted_bam}
        samtools index -@ {threads} {output.sorted_bam}
        echo "END indexing"
        echo "END indexing" 1>&2

        samtools idxstats {output.sorted_bam} > {output.idxstat}
        echo "END idxstats"
        echo "END idxstats" 1>&2
        """

rule convert_bam:
    input:
        bam="analysis/bowtie2/{sample}.bam"
    output:
        sorted_bed="analysis/bed_files/{sample}.bed",
        gzip_bed="analysis/bed_files/{sample}.bed.gz"
    benchmark:
        "benchmarks/bed_files/{sample}.txt"
    envmodules:
        config['modules']['bedtools'],
        config['modules']['htslib']
    threads: 4
    resources:
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        bedtools bamtobed -bedpe -i {input.bam} | \
        cut -f1,2,6,7 | \
        awk -F'\t' '$1 != "." || $2 != -1 || $3 != -1' | \
        sort -k1,1 -k2,2n -k3,3n > {output.sorted_bed}

        bgzip -c {output.sorted_bed} > {output.gzip_bed}
        """

rule generate_bw:
    input:
        gzip_bed="analysis/bed_files/{sample}.bed.gz"
    output:
        bw_file="analysis/bigwig_files/{sample}.bw",
        bedgraph="analysis/bigwig_files/{sample}.bedgraph"
    benchmark:
        "benchmarks/bigwig_files/{sample}.txt"
    envmodules:
        config['modules']['bedtools'],
        config['modules']['ucsc']
    params:
        chr_len=config['ref']['fai']
    threads: 4
    resources:
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        coverage=$(zcat {input.gzip_bed} | awk '{{s+=($3-$2)}} END {{print s}}')
        echo "Coverage: $coverage"

        scaling_factor=$(awk -v c="$coverage" 'BEGIN {{print (1/c)*10^10}}')
        echo "Scaling Factor: $scaling_factor"

        zcat {input.gzip_bed} | bedtools genomecov -bg -i stdin -g {params.chr_len} -scale $scaling_factor > {output.bedgraph}
        echo "Bedgraph file created: {output.bedgraph}"

        bedGraphToBigWig {output.bedgraph} {params.chr_len} {output.bw_file}
        echo "BigWig file created: {output.bw_file}"

        """
