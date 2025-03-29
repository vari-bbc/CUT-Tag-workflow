import pandas as pd
import numpy as np
import os
import re
import itertools
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("7.25.0")

##### Editable variables #####

configfile_path = "bin/config.yaml"
configfile: configfile_path

bt2_index = config['ref']['index']

##### load config and sample sheets #####
samplesheet="bin/samples.tsv"
units = pd.read_table(samplesheet, dtype={"sample" : str, "sample_group" : str })

validate(units, "schema/samplesheet.yaml")

contrasts_file="bin/contrasts.tsv"
contrasts = pd.read_table(contrasts_file, dtype={"group1" : str, "group2" : str , "enriched_factor" : str})
validate(contrasts, "schema/contrasts.yaml")
contrasts["enriched_factor"] = contrasts["enriched_factor"].fillna("peaks")

samples = units[["sample","control","sample_group","enriched_factor"]].drop_duplicates()
if not samples['sample'].is_unique:
    raise Exception('A sample has more than one combination of control, sample_group, enriched_factor.')

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

bigwig_norms = ['baseCov']
bigwig_norms = bigwig_norms + config['addtnl_bigwig_norms'] if isinstance(config['addtnl_bigwig_norms'], list) else bigwig_norms

rule all:
    input:
        "analysis/multiqc/multiqc_report.html",
        expand("analysis/bigwig_files/{norm_method}/{sample.sample}.bw", sample=samples.itertuples(), norm_method = [k for k in bigwig_norms if not 'csaw' in k]),
        expand("analysis/bigwig_files/{norm_method}/{sample.sample}.bw", sample=samples_no_controls.itertuples(), norm_method = [k for k in bigwig_norms if 'csaw' in k]),
        expand("analysis/macs3_narrow/{sample.sample}_summits.bed", sample=samples.itertuples()),
        expand("analysis/macs3_broad/{sample.sample}_peaks.broadPeak", sample=samples.itertuples()) if config['macs3']['run_broad'] else [],
        #expand("analysis/{peak_type}/merged/{merge_id}.bed", peak_type = peak_types, merge_id = sample_groups + ['all']),
        #"analysis/csaw_count/peaks/global_filt.rds",
        expand("analysis/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/results_se.rds", enriched_factor=pd.unique(samples_no_controls['enriched_factor']), peak_type = config['csaw']['peak_type'], norm_type = config['csaw']['norm_type']),
        expand("analysis/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/csaw_summary.html", enriched_factor=pd.unique(samples_no_controls['enriched_factor']), peak_type = config['csaw']['peak_type'], norm_type = config['csaw']['norm_type'])

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
        unsorted_bam="analysis/bowtie2/{sample}.bam",
        outbai="analysis/bowtie2/{sample}.sorted.bam.bai",
    benchmark:
        "benchmarks/bowtie2/{sample}.txt"
    params:
        bt2_index=bt2_index,
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
        bowtie2 -p {threads} --very-sensitive-local --soft-clipped-unmapped-tlen --no-mixed --no-discordant --dovetail --phred33 \
        -I 10 -X 1000 -x {params.bt2_index} -1 {input.R1_trimmed} -2 {input.R2_trimmed} | \
        samblaster | \
        samtools view -bS \
        -@ {threads} \
        -O "BAM" \
        -o {output.unsorted_bam} 
        
        samtools sort -m 6G -O "bam" -@ {threads} -o {output.sorted_bam} {output.unsorted_bam}
        samtools index -@ {threads} {output.sorted_bam}

        """


rule get_keep_chrom_names:
    input:
        ref_fasta=config["ref"]["sequence"],
    output:
        "analysis/misc/keep_chroms.txt"
    benchmark:
        "benchmarks/get_keep_chrom_names/bench.txt"
    params:
        keep_std_chroms=config["keep_std_chroms"],
        rm_chroms="foobar" if config["rm_chroms"] == "" else config["rm_chroms"],
    threads: 1
    resources:
        mem_gb=20,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
        config['modules']['R']
    script:
        "bin/scripts/get_keep_chrom_names.R"

def get_keep_regions (wildcards):
    cmd = ""
    if config['ref']['blacklist'] == "":
        cmd = "perl -F/\t/ -lane qq|print $F[0]\t0\t$F[1]\n|" + config['ref']['fai']
    else:
        cmd = "bedtools sort -i " + config['ref']['blacklist'] + " -g " + config['ref']['fai'] + " | bedtools complement -i 'stdin' -g " + config['ref']['fai']

    return cmd

rule make_keep_regions_bed:
    """
    Define regions to keep alignments for. First, find the complement of blacklist regions (if provided) then keep only standard chromosomes (removing certain chromosomes if requested).
    """
    input:
        keep_chroms = "analysis/misc/keep_chroms.txt" if config['keep_chroms_file'] == "" else config['keep_chroms_file'],
    output:
        keep_regions_temp = "analysis/misc/keep_regions.bed.temp",
        keep_regions = "analysis/misc/keep_regions.bed"
    params:
        keep_regions = get_keep_regions,
    benchmark:
        "benchmarks/make_keep_regions_bed/bench.txt"
    envmodules:
        config['modules']['bedtools'],
    threads: 8
    resources:
        mem_gb = 96,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        {params.keep_regions} > {output.keep_regions_temp}

        for chr in $(cat {input.keep_chroms})
        do
            grep -P "^$chr\t" {output.keep_regions_temp} >> {output.keep_regions} 
        done

        """

rule filter_bams:
    """
    Filter BAMs. After filters, fixmate will be run to update SAM flags for pairs where one mate is removed and then filtered for properly paired alignments
    """
    input:
        bam = "analysis/{align_dirname}/{bam_name}.sorted.bam",
        keep_regions = "analysis/misc/keep_regions.bed"
    output:
        sorted_bam="analysis/{align_dirname}_filt/{bam_name}.sorted.bam",
        sorted_bai="analysis/{align_dirname}_filt/{bam_name}.sorted.bam.bai",
        bam="analysis/{align_dirname}_filt/{bam_name}.bam",
    params:
        view_mapq="" if config['samtools_mapq'] == "" else "-q {mapq}".format(mapq=config['samtools_mapq']),
        view_keep="" if config['samtools_keep_flags'] == "" else "-f {flag}".format(flag=config['samtools_keep_flags']),
        view_omit="" if config['samtools_omit_flags'] == "" else "-F {flag}".format(flag=config['samtools_omit_flags']),
    benchmark:
        "benchmarks/{align_dirname}_filt/{bam_name}.txt"
    envmodules:
        config['modules']['samtools']
    threads: 8
    resources:
        mem_gb = 96,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        samtools view -@ {threads} -u --region-file {input.keep_regions} {params.view_mapq} {params.view_keep} {params.view_omit} {input.bam} | \
        samtools sort -@ {threads} -n -u - | \
        samtools fixmate -@ {threads} -u - - | \
        samtools view -@ {threads} -f 2 -o {output.bam} -

        samtools sort -@ {threads} -o {output.sorted_bam} {output.bam}
        samtools index -@ {threads} {output.sorted_bam}

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

rule qualimap:
    """
    Run qualimap on unfiltered alignments.
    """
    input:
        "analysis/{align_dirname}/{sample}.sorted.bam",
    output:
        done=touch("analysis/{align_dirname}/qualimap/{sample}/done")
    benchmark:
        "benchmarks/{align_dirname}/qualimap/{sample}.txt"
    envmodules:
        config['modules']['qualimap']
    params:
        outdir=lambda wildcards,output: os.path.dirname(output.done)
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards)
    threads: 8
    shell:
        """
        qualimap bamqc -bam {input} --java-mem-size={resources.mem_gb}G --paint-chromosome-limits -outdir {params.outdir} -nt {threads}

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
    """
    Assumes non-discordant alignments as it outputs a BED file that represents the alignment of the fragment not the individual reads.
    """
    input:
        bam="analysis/bowtie2_filt/{sample}.bam"
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

rule base_cov_scale_factors:
    """
    Generate scaling factor for generate_bw rule to normalize by total base coverage.
    """
    input:
        gzip_bed="analysis/bed_files/{sample}.bed.gz"
    output:
        "analysis/bigwig_norm_factors/base_cov_scale_factors/{sample}.tsv"
    benchmark:
        "benchmarks/bigwig_norm_factors/base_cov_scale_factors/{sample}.txt"
    params:
    threads: 4
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
    shell:
        """
        coverage=$(zcat {input.gzip_bed} | awk '{{s+=($3-$2)}} END {{print s}}')
        echo "Coverage: $coverage"

        scaling_factor=$(awk -v c="$coverage" 'BEGIN {{print (1/c)*10^10}}')
        echo "Scaling Factor: $scaling_factor"

        echo "$scaling_factor" > {output}
        """

rule concat_base_cov_scale_factors:
    """
    Concatenate the total base coverage normalization factors.
    """
    input:
        expand("analysis/bigwig_norm_factors/base_cov_scale_factors/{sample}.tsv", sample = samples['sample'])
    output:
        "analysis/bigwig_norm_factors/base_cov_scale_factors.tsv"
    benchmark:
        "benchmarks/bigwig_norm_factors/concat_base_cov_scale_factors.tsv"
    params:
    threads: 1
    resources:
        mem_gb=32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
    shell:
        """
        printf "sample\\tfinal.factors\\n" > {output}
        for tsv in {input}
        do
            printf "$(basename $tsv | perl -npe 's:.tsv::')\\t$(head -n1 $tsv)\\n" >> {output}
        done
        """

def get_scale_factor_file(wildcards):
    curr_enriched = ""
    if (wildcards.sample in samples_no_controls['sample'].values):
        curr_enriched = samples_no_controls[samples_no_controls['sample']==wildcards.sample]['enriched_factor'].values[0]
    if (wildcards.scale_method == "csaw_bkgd"):
        return "analysis/bigwig_norm_factors/{enriched_factor}_csaw_bkgd.tsv".format(enriched_factor=curr_enriched)
    elif (wildcards.scale_method == "csaw_hiAbund"):
        return "analysis/bigwig_norm_factors/{enriched_factor}_csaw_hiAbund.tsv".format(enriched_factor=curr_enriched)
    elif (wildcards.scale_method == "baseCov"):
        return "analysis/bigwig_norm_factors/base_cov_scale_factors.tsv"
    else:
        raise Exception("Could not determine scale factors file name.")

def get_bigwig_norm_factor(wildcards, input):
    df = pd.read_table(input.scale_factor)
    scalefactor = df[df['sample']==wildcards.sample]['final.factors'].values[0]
    return str(scalefactor)

rule generate_bw:
    input:
        gzip_bed="analysis/bed_files/{sample}.bed.gz",
        scale_factor=get_scale_factor_file
    output:
        bw_file="analysis/bigwig_files/{scale_method}/{sample}.bw",
        bedgraph=temp("analysis/bigwig_files/{scale_method}/{sample}.bedgraph")
    benchmark:
        "benchmarks/bigwig_files/{scale_method}/{sample}.txt"
    envmodules:
        config['modules']['bedtools'],
        config['modules']['ucsc']
    params:
        chr_len=config['ref']['fai'],
        scale_factor=get_bigwig_norm_factor,
    threads: 4
    resources:
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        scaling_factor="{params.scale_factor}"
        echo "Scaling Factor: $scaling_factor"

        zcat {input.gzip_bed} | bedtools genomecov -bg -i stdin -g {params.chr_len} -scale $scaling_factor > {output.bedgraph}
        echo "Bedgraph file created: {output.bedgraph}"

        bedGraphToBigWig {output.bedgraph} {params.chr_len} {output.bw_file}
        echo "BigWig file created: {output.bw_file}"

        """


def get_macs3_bams(wildcards):
    macs3_bams = { 'trt': "analysis/bowtie2_filt/{sample}.sorted.bam".format(sample=wildcards.sample) }
    
    control = samples[samples['sample'] == wildcards.sample]['control'].values[0]
    
    if (not pd.isnull(control)):
        macs3_bams['control'] = expand("analysis/bowtie2_filt/{sample}.sorted.bam", sample = control.split(','))
    
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
    df = samples_no_controls
    if (wildcards.merge_id in df.sample_group.values):
        curr_samples = df[(df['sample_group'] == wildcards.merge_id) & (df['enriched_factor'] == wildcards.enriched_factor)]['sample'].values
    elif (wildcards.merge_id == "all"):
        curr_samples = df[df['enriched_factor'] == wildcards.enriched_factor]['sample'].values
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
        "analysis/{peak_type}/merged_{enriched_factor}/{merge_id}.bed"
    benchmark:
        "benchmarks/{peak_type}/merged_{enriched_factor}/{merge_id}.txt"
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

rule frags_over_peaks:
    input:
        frags = "analysis/bed_files/{sample}.bed.gz",
        peaks = lambda wildcards: "analysis/{{peak_type}}/merged_{enriched_factor}/{{merge_id}}.bed".format(enriched_factor = samples_no_controls[samples_no_controls['sample'] == wildcards.sample]['enriched_factor'].values[0])
    output:
        bed = "analysis/frags_over_peaks/{peak_type}/{merge_id}/{sample}.bed",
        frag_count = "analysis/frags_over_peaks/{peak_type}/{merge_id}/{sample}.txt"
    benchmark:
        "benchmarks/frags_over_peaks/{peak_type}/{merge_id}/{sample}.txt"
    params:
    envmodules:
        config['modules']['bedtools']
    threads: 1
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        tot_frags=`zcat {input.frags} | wc -l`
        overlap_frags=`bedtools intersect -u -a {input.frags} -b {input.peaks} | tee {output.bed} | wc -l`
        frip=`echo "scale=5;$overlap_frags / $tot_frags" | bc`

        printf "$tot_frags\\n$overlap_frags\\n$frip\\n" > {output.frag_count}
        """

def get_frag_count_files (wildcards):
    out_files = expand("analysis/frags_over_peaks/{peak_type}/{sample.sample_group}/{sample.sample}.txt", peak_type = config['frip']['peak_type'], sample = samples_no_controls.itertuples())
    out_files = out_files + expand("analysis/frags_over_peaks/{peak_type}/all/{sample.sample}.txt", peak_type = config['frip']['peak_type'], sample = samples_no_controls.itertuples())
    return out_files

rule frip:
    input:
        frag_counts=get_frag_count_files,
        samplesheet=samplesheet
    output:
        all_merge="analysis/frags_over_peaks/all_frip.csv",
        group_merge="analysis/frags_over_peaks/group_frip.csv"
    benchmark:
        "benchmarks/frags_over_peaks/frip.txt"
    params:
        peak_type=config['frip']['peak_type'],
        noncontrol_samples=','.join(samples_no_controls['sample'].values)
    envmodules:
        config['modules']['R']
    threads: 1
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    script:
        "bin/scripts/calc_frip.R"


rule csaw_count:
    input:
        samplesheet=samplesheet,
        bams = lambda wildcards: expand("analysis/bowtie2_filt/{sample}.sorted.bam", sample=samples_no_controls[samples_no_controls['enriched_factor'] == wildcards.enriched_factor]['sample']),
    output:
        binned="analysis/csaw_count/{enriched_factor}/binned.rds",
        small_wins="analysis/csaw_count/{enriched_factor}/small_wins.rds",
        filt_small_wins="analysis/csaw_count/{enriched_factor}/filt_small_wins.rds",
        global_filt="analysis/csaw_count/{enriched_factor}/global_filt.rds",
        bkgrd_scale="analysis/csaw_count/{enriched_factor}/bkgrd_norm_factors.tsv",
        hiAbund_scale="analysis/csaw_count/{enriched_factor}/hiAbund_norm_factors.tsv",
        bkgrd_scale_link="analysis/bigwig_norm_factors/{enriched_factor}_csaw_bkgd.tsv",
        hiAbund_scale_link="analysis/bigwig_norm_factors/{enriched_factor}_csaw_hiAbund.tsv",
        figs_dir=directory("analysis/csaw_count/{enriched_factor}/figures")
    benchmark:
        "benchmarks/csaw_count/{enriched_factor}.txt"
    params:
        window_width=config['csaw']['win_width'],
        samp_names=lambda wildcards, input: [os.path.basename(x).replace(".sorted.bam","") for x in input.bams],
    threads: 16
    resources:
        mem_gb=396,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
        config['modules']['R']
    script:
        "bin/scripts/csaw_count.R"

rule csaw_diff:
    input:
        filt_small_wins="analysis/csaw_count/{enriched_factor}/filt_small_wins.rds",
        bin_counts="analysis/csaw_count/{enriched_factor}/binned.rds",
        merged_peaks="analysis/{peak_type}/merged_{enriched_factor}/all.bed",
        contrasts="bin/contrasts.tsv"
    output:
        edgeR_rds="analysis/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/edgeR_objs.rds",
        results_rds="analysis/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/results.rds",
        results_se_rds="analysis/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/results_se.rds",
        peaks_res_rds="analysis/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/peaks_res.rds",
        combined_rds="analysis/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/combined.rds",
        peaks_anno_rds="analysis/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/peaks_anno.rds",
        figs_dir=directory("analysis/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/figures")
    benchmark:
        "benchmarks/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}.txt"
    params:
        orgdb=config['chipseeker_params']['orgdb'],
        txdb=config['chipseeker_params']['txdb'],
        promo_start=config['chipseeker_params']['promoter_start'],
        promo_end=config['chipseeker_params']['promoter_end']
    threads: 16
    resources:
        mem_gb=396,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
        config['modules']['R']
    script:
        "bin/scripts/csaw_diff.R"

rule csaw_summary:
    input:
        rmd="bin/scripts/csaw_summary.Rmd",
        edgeR_rds="analysis/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/edgeR_objs.rds",
        peaks_res_rds="analysis/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/peaks_res.rds",
        combined_rds="analysis/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/combined.rds",
    output:
        rmd="analysis/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/csaw_summary.Rmd",
        edgeR_rds="analysis/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/edgeR_objs.rds",
        peaks_res_rds="analysis/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/peaks_res.rds",
        combined_rds="analysis/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/combined.rds",
        out_res="analysis/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/DB_results.xlsx",
        html_report="analysis/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/csaw_summary.html",
        figdir=directory("analysis/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/individual_figures")
    benchmark:
        "benchmarks/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}.txt"
    params:
        out_res=lambda wildcards, output: os.path.basename(output.out_res),
        fdr_th = "0.05",
        wd=lambda wildcards, output: os.path.dirname(output.html_report),
        rmd = lambda wildcards, output: os.path.basename(output.rmd),
        edgeR_rds = lambda wildcards, output: os.path.basename(output.edgeR_rds),
        peaks_res_rds = lambda wildcards, output: os.path.basename(output.peaks_res_rds),
        combined_rds = lambda wildcards, output: os.path.basename(output.combined_rds),
        figdir = lambda wildcards, output: os.path.basename(output.figdir)
    threads: 4
    resources:
        mem_gb=64,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
        config['modules']['R'],
        config['modules']['pandoc']
    shell:
        """
        cp {input.rmd} {output.rmd}
        cp {input.edgeR_rds} {output.edgeR_rds}
        cp {input.peaks_res_rds} {output.peaks_res_rds}
        cp {input.combined_rds} {output.combined_rds}

        cd {params.wd}

        Rscript --vanilla -e "rmarkdown::render('{params.rmd}', params = list(out_res = '{params.out_res}', fdr_th = '{params.fdr_th}', edgeR_rds = '{params.edgeR_rds}', peaks_res_rds = '{params.peaks_res_rds}', combined_rds = '{params.combined_rds}', figdir = '{params.figdir}'))"
        """

rule multiqc:
    input:
        expand("analysis/fastqc/{sample.sample}_R1_fastqc.html", sample=samples.itertuples()),
        expand("analysis/fastqc/{sample.sample}_R2_fastqc.html", sample=samples.itertuples()),
        expand("analysis/trim_galore/{sample.sample}_R{read}_val_{read}_fastqc.html", sample=samples.itertuples(), read=["1","2"]),
        expand("analysis/{align_dir}/idxstats/{sample.sample}.idxstats", sample=samples.itertuples(), align_dir=['bowtie2','bowtie2_filt']),
        expand("analysis/{align_dir}/flagstat/{sample.sample}.flagstat", sample=samples.itertuples(), align_dir=['bowtie2','bowtie2_filt']),
        expand("analysis/{align_dir}/CollectAlignmentSummaryMetrics/{sample.sample}.aln_metrics.txt", sample=samples.itertuples(), align_dir=['bowtie2','bowtie2_filt']),
        expand("analysis/{align_dir}/CollectInsertSizeMetrics/{sample.sample}.insert_size_metrics.txt", sample=samples.itertuples(), align_dir=['bowtie2','bowtie2_filt']),
        expand("analysis/{align_dir}/qualimap/{sample.sample}/done", sample=samples.itertuples(), align_dir=['bowtie2','bowtie2_filt']),
        expand("analysis/macs3_narrow/{sample.sample}_peaks.xls", sample=samples.itertuples()),
        expand("analysis/frags_over_peaks/{peak_type}_frip.csv", peak_type = ['all','group'])
    output:
        mqc_vers="analysis/multiqc/all_mqc_versions.yaml",
        mqc="analysis/multiqc/multiqc_report.html",
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    params:
        config_file=configfile_path,
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
        grep -P -A999 "^modules:"  {params.config_file} | tail -n+2 | perl -npe 's:bbc2/[^/]+/[^\-]+\-::'  > {output.mqc_vers}

        multiqc \
        --force \
        --config bin/multiqc_config2.yaml \
        --outdir {params.workdir} \
        --filename {params.outfile} \
        {params.dirs} analysis/multiqc/all_mqc_versions.yaml

        """
