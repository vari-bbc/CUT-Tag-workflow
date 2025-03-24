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

# Filter for sample rows that are not controls
controls_list = list(itertools.chain.from_iterable( [x.split(',') for x in samples['control'].values if not pd.isnull(x)] ))
samples_no_controls = samples[-samples['sample'].isin(controls_list)].copy()
samples_no_controls["enriched_factor"] = samples_no_controls["enriched_factor"].fillna("peaks")
sample_groups = pd.unique(samples_no_controls['sample_group']).tolist()

snakemake_dir = os.getcwd() + "/"

# make a tmp directory for analyses
tmp_dir = os.path.join(snakemake_dir, "tmp")
if not os.path.exists(tmp_dir):
    os.mkdir(tmp_dir)

peak_types = ['macs3_narrow']
peak_types = peak_types + ["macs3_broad"] if config['macs3']['run_broad'] else peak_types

rule all:
    input:
        "analysis/multiqc/multiqc_report.html",
        expand("analysis/bigwig_files/{sample.sample}.bw", sample=samples.itertuples()),
        expand("analysis/macs3_narrow/{sample.sample}_summits.bed", sample=samples.itertuples()),
        expand("analysis/macs3_broad/{sample.sample}_peaks.broadPeak", sample=samples.itertuples()) if config['macs3']['run_broad'] else [],
        expand("analysis/{peak_type}/merged/{merge_id}.bed", peak_type = peak_types, merge_id = sample_groups + ['all']),

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

rule trim_galore_PE:
    """
    Run trim_galore on paired-end reads.
    """
    input:
        expand("analysis/renamed_data/{{sample}}_R{read}.fastq.gz", read=["1","2"])
    output:
        temp(expand("analysis/trim_galore/{{sample}}_R{ext}", ext=["1_val_1.fq.gz","2_val_2.fq.gz"])),
        expand("analysis/trim_galore/{{sample}}_R1{ext}", ext=[".fastq.gz_trimming_report.txt","_val_1_fastqc.html"]),
        expand("analysis/trim_galore/{{sample}}_R2{ext}", ext=[".fastq.gz_trimming_report.txt","_val_2_fastqc.html"]),
    params:
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        trim_galore --paired {input} --output_dir analysis/trim_galore/ --cores {threads} --fastqc
        """


rule bowtie2:
    input:
        R1_trimmed="analysis/trim_galore/{sample}_R1_val_1.fq.gz",
        R2_trimmed="analysis/trim_galore/{sample}_R2_val_2.fq.gz"   
    output:
        sorted_bam="analysis/bowtie2/{sample}.sorted.bam",
        unsorted_bam=temp("analysis/bowtie2/{sample}.bam"),
        outbai="analysis/bowtie2/{sample}.sorted.bam.bai",
    benchmark:
        "benchmarks/bowtie2/{sample}.txt"
    params:
        bt2_index=bt2_index,
        samblaster_params=lambda wildcards: "--addMateTags" if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "--ignoreUnmated",
        dedup_params="--removeDups" if config['dedup'] else "",
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
        -I 10 -X 1000 -x {params.bt2_index} -1 {input.R1_trimmed} -2 {input.R2_trimmed} | \
        samblaster {params.dedup_params} | \
        samtools view -bS \
        -@ {threads} \
        -O "BAM" \
        -o {output.unsorted_bam} 
        
        samtools sort {output.unsorted_bam} -m 6G -O "bam" -o {output.sorted_bam}
        samtools index -@ {threads} {output.sorted_bam}
        echo "END indexing"
        echo "END indexing" 1>&2

        """


rule idxstats:
    """
    Run samtools idxstats.
    """
    input:
        "analysis/{align_dirname}/{bam_name}.sorted.bam"
    output:
        "analysis/{align_dirname}/idxstats/{bam_name}.idxstats"
    params:
    benchmark:
        "benchmarks/{align_dirname}/idxstats/{bam_name}.txt"
    envmodules:
        config['modules']['samtools']
    threads: 8
    resources:
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        samtools idxstats -@ {threads} {input} > {output}
        """

rule flagstat:
    """
    Run samtools flagstat.
    """
    input:
        "analysis/{align_dirname}/{bam_name}.sorted.bam"
    output:
        "analysis/{align_dirname}/flagstat/{bam_name}.flagstat"
    params:
    benchmark:
        "benchmarks/{align_dirname}/flagstat/{bam_name}.txt"
    envmodules:
        config['modules']['samtools']
    threads: 8
    resources:
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        samtools flagstat -@ {threads} {input} > {output}
        """

rule CollectAlignmentSummaryMetrics:
    """
    Run Picard CollectAlignmentSummaryMetrics.
    """
    input:
        "analysis/{align_dirname}/{bam_name}.sorted.bam"
    output:
        out="analysis/{align_dirname}/CollectAlignmentSummaryMetrics/{bam_name}.aln_metrics.txt",
    params:
        temp="./analysis/{align_dirname}/CollectAlignmentSummaryMetrics/",
        reffasta=config["ref"]["sequence"]
    benchmark:
        "benchmarks/{align_dirname}/CollectAlignmentSummaryMetrics/{bam_name}.txt"
    envmodules:
        config['modules']['picard']
    threads: 4
    resources:
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir={params.temp} -jar $PICARD CollectAlignmentSummaryMetrics I={input} O={output.out} R={params.reffasta}
        """

rule CollectInsertSizeMetrics:
    """
    Run Picard CollectInsertSizeMetrics.
    """
    input:
        "analysis/{align_dirname}/{bam_name}.sorted.bam"
    output:
        out="analysis/{align_dirname}/CollectInsertSizeMetrics/{bam_name}.insert_size_metrics.txt",
        hist="analysis/{align_dirname}/CollectInsertSizeMetrics/{bam_name}.insert_size_histogram.pdf"
    params:
        temp="./analysis/{align_dirname}/CollectInsertSizeMetrics/"
    benchmark:
        "benchmarks/{align_dirname}/CollectInsertSizeMetrics/{bam_name}.txt"
    envmodules:
        config['modules']['picard']
    threads: 4
    resources:
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir={params.temp} -jar $PICARD CollectInsertSizeMetrics I={input} O={output.out} H={output.hist}
        """


rule convert_bam:
    input:
        bam="analysis/bowtie2/{sample}.bam"
    output:
        sorted_bed=temp("analysis/bed_files/{sample}.bed"),
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
        bedgraph=temp("analysis/bigwig_files/{sample}.bedgraph")
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


def get_macs3_bams(wildcards):
    macs3_bams = { 'trt': "analysis/bowtie2/{sample}.sorted.bam".format(sample=wildcards.sample) }
    
    control = samples[samples['sample'] == wildcards.sample]['control'].values[0]
    
    if (not pd.isnull(control)):
        macs3_bams['control'] = expand("analysis/bowtie2/{sample}.sorted.bam", sample = control.split(','))
    
    return macs3_bams


rule macs3_narrow:
    shadow: "shallow" # macs3 makes a tmp_b_c.txt file in the working directory; not sure how that is handled when multiple instances of macs3 are run simultaneously. To be safe, will run using a shadow rule so that each instance can create its own tmp_b_c.txt
    input:
        unpack(get_macs3_bams)
    output:
        multiext("analysis/macs3_narrow/{sample}_peaks", ".xls", ".narrowPeak"),
        "analysis/macs3_narrow/{sample}_summits.bed",
    benchmark:
        "benchmarks/macs3_narrow/{sample}.txt"
    params:
        control_param=lambda wildcards, input: "-c " + ' '.join(input.control) if 'control' in input.keys() else '',
        #peak_type_param=lambda wildcards: "--broad" if wildcards.peak_type == "broad" else "",
        species=config['macs3']['species'],
        q_cutoff="0.01",
        name="{sample}",
        outdir="analysis/macs3_narrow/",
        tmpdir=tmp_dir,
    envmodules:
        config['modules']['macs3']
    threads: 1
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """

        macs3 callpeak \
        -t {input.trt} \
        {params.control_param} \
        -f BAMPE \
        -g {params.species} \
        -n {params.name} \
        -q {params.q_cutoff} \
        --keep-dup all \
        --outdir {params.outdir} \
        --tempdir {params.tmpdir}

        """

rule macs3_broad:
    shadow: "shallow" # macs3 makes a tmp_b_c.txt file in the working directory; not sure how that is handled when multiple instances of macs3 are run simultaneously. To be safe, will run using a shadow rule so that each instance can create its own tmp_b_c.txt
    input:
        unpack(get_macs3_bams)
    output:
        multiext("analysis/macs3_broad/{sample}_peaks", ".xls", ".broadPeak", ".gappedPeak"),
    benchmark:
        "benchmarks/macs3_broad/{sample}.txt"
    params:
        control_param=lambda wildcards, input: "-c " + ' '.join(input.control) if 'control' in input.keys() else '',
        #peak_type_param=lambda wildcards: "--broad" if wildcards.peak_type == "broad" else "",
        species=config['macs3']['species'],
        q_cutoff="0.01",
        name="{sample}",
        outdir="analysis/macs3_broad/",
        tmpdir=tmp_dir,
    envmodules:
        config['modules']['macs3']
    threads: 1
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """

        macs3 callpeak \
        --broad \
        -t {input.trt} \
        {params.control_param} \
        -f BAMPE \
        -g {params.species} \
        -n {params.name} \
        -q {params.q_cutoff} \
        --keep-dup all \
        --outdir {params.outdir} \
        --tempdir {params.tmpdir}

        """

def get_beds_for_merging (wildcards):
    curr_samples = []
    if (wildcards.merge_id in samples_no_controls.sample_group.values):
        curr_samples = samples_no_controls[samples_no_controls['sample_group'] == wildcards.merge_id]['sample'].values
    elif (wildcards.merge_id == "all"):
        curr_samples = samples_no_controls['sample'].values
    else:
        raise Exception("Invalid merge_id for peak merging.")

    in_beds = []
    if (wildcards.peak_type == "macs3_narrow"):
        in_beds = expand("analysis/macs3_narrow/{sample}_peaks.narrowPeak", sample = curr_samples)
    elif (wildcards.peak_type == "macs3_broad"):
        in_beds = expand("analysis/macs3_broad/{sample}_peaks.broadPeak", sample = curr_samples)
    else:
        raise Exception("Invalid peak_type for peak merging.")
    
    return in_beds

rule merge_peaks:
    input:
        get_beds_for_merging
    output:
        "analysis/{peak_type}/merged/{merge_id}.bed"
    benchmark:
        "benchmarks/{peak_type}/merged/{merge_id}.txt"
    params:
    envmodules:
        config['modules']['bedtools']
    threads: 1
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """

        cat {input} | sort -k 1,1 -k2,2n | bedtools merge > {output}
        """

rule multiqc:
    input:
        expand("analysis/fastqc/{sample.sample}_R1_fastqc.html", sample=samples.itertuples()),
        expand("analysis/fastqc/{sample.sample}_R2_fastqc.html", sample=samples.itertuples()),
        expand("analysis/trim_galore/{sample.sample}_R{read}_val_{read}_fastqc.html", sample=samples.itertuples(), read=["1","2"]),
        expand("analysis/bowtie2/idxstats/{sample.sample}.idxstats", sample=samples.itertuples()),
        expand("analysis/bowtie2/flagstat/{sample.sample}.flagstat", sample=samples.itertuples()),
        expand("analysis/bowtie2/CollectAlignmentSummaryMetrics/{sample.sample}.aln_metrics.txt", sample=samples.itertuples()),
        expand("analysis/bowtie2/CollectInsertSizeMetrics/{sample.sample}.insert_size_metrics.txt", sample=samples.itertuples()),
        expand("analysis/macs3_narrow/{sample.sample}_peaks.xls", sample=samples.itertuples()),
    output:
        "analysis/multiqc/multiqc_report.html",
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    params:
        workdir="analysis/multiqc/",
        dirs=lambda wildcards,input: " ".join(pd.unique([os.path.dirname(x) for x in input])),
        outfile="multiqc_report"
    envmodules:
        config['modules']['multiqc']
    threads: 4
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        multiqc \
        --force \
        --outdir {params.workdir} \
        --filename {params.outfile} \
        {params.dirs}

        """
