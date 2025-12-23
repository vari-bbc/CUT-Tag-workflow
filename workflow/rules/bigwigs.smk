rule convert_bam:
    """
    Assumes non-discordant alignments as it outputs a BED file that represents the alignment of the fragment not the individual reads.
    """
    input:
        bam="analysis/bowtie2_filt_endogenous/{sample}.bam"
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
        return "analysis/bigwig_norm_factors/{enriched_factor}_csaw_win{width}_bkgd_endogenous.tsv".format(enriched_factor=curr_enriched, width=csaw_win_sizes[0])
    elif (wildcards.scale_method == "csaw_hiAbund"):
        return "analysis/bigwig_norm_factors/{enriched_factor}_csaw_win{width}_hiAbund_endogenous.tsv".format(enriched_factor=curr_enriched, width=csaw_win_sizes[0])
    elif (wildcards.scale_method == "baseCov"):
        return "analysis/bigwig_norm_factors/base_cov_scale_factors.tsv"
    elif (wildcards.scale_method == "csaw.bkgd_spikein"):
        return "analysis/bigwig_norm_factors/{enriched_factor}_csaw_win{width}_bkgd_spikein.tsv".format(enriched_factor=curr_enriched, width=csaw_win_sizes[0])
    elif (wildcards.scale_method == "csaw.hiAbund_spikein"):
        return "analysis/bigwig_norm_factors/{enriched_factor}_csaw_win{width}_hiAbund_spikein.tsv".format(enriched_factor=curr_enriched, width=csaw_win_sizes[0])
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


rule avg_bigwigs:
    input:
        bw=lambda wildcards: expand("analysis/bigwig_files/{scale_method}/{sample}.bw", sample=samples[samples['sample_group']==wildcards.group]['sample'].values, scale_method=wildcards.scale_method)
    output:
        bw="analysis/avg_bigwigs/{scale_method}/{group}.bw"
    params:
    benchmark:
        "benchmarks/avg_bigwigs/{scale_method}/{group}.txt"
    envmodules:
        config['modules']['deeptools']
    threads: 8
    resources:
        nodes = 1,
        mem_gb = 72,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        bigwigAverage -b {input.bw} --binSize 1 -p {threads} -o {output.bw} -of "bigwig"
        """