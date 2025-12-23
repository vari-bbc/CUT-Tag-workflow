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

rule make_mqc_sample_filters:
    input:
        samplesheet=samplesheet
    output:
        sample_filts="analysis/misc/mqc_sample_filters.tsv"
    benchmark:
        "benchmarks/make_mqc_sample_filters/bench.txt"
    params:
        noncontrol_samples=','.join(samples_no_controls['sample'].values)
    envmodules:
        config['modules']['R']
    threads: 1
    resources:
        mem_gb=8,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    script:
        "bin/scripts/make_mqc_sample_filters.R"

rule make_mqc_config:
    input:
        "bin/multiqc_config.yaml",
    output:
        "analysis/misc/mqc_config.yaml"
    benchmark:
        "benchmarks/make_mqc_config/bench.txt"
    params:
        wd=snakemake_dir
    envmodules:
    threads: 1
    resources:
        mem_gb=8,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        cp {input} {output}
        perl -i -lnpe 's:images/VAI_2_Line_White.png:{params.wd}images/VAI_2_Line_White.png:' {output}
        """

rule fastq_screen:
    """
    Run fastq_screen to detect any contamination from other species or excessive rRNA.
    """
    input:
        "analysis/renamed_data/{fq_pref}.fastq.gz"
    output:
        html = "analysis/fastq_screen/{fq_pref}_screen.html",
        txt = "analysis/fastq_screen/{fq_pref}_screen.txt",
    params:
    benchmark:
        "benchmarks/fastq_screen/{fq_pref}.txt"
    envmodules:
        config['modules']['fastq_screen']
    threads: 8
    resources:
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        fastq_screen --threads {threads} --outdir analysis/fastq_screen/ {input}
        """


rule multiqc:
    input:
        mqc_config="analysis/misc/mqc_config.yaml",
        sample_filters="analysis/misc/mqc_sample_filters.tsv",
        mqc_dirs=expand("analysis/fastqc/{sample.sample}_R1_fastqc.html", sample=samples.itertuples()) +
            expand("analysis/fastqc/{sample.sample}_R2_fastqc.html", sample=samples.itertuples()) +
            expand("analysis/fastq_screen/{sample.sample}_R1_screen.html", sample=samples.itertuples()) +
            expand("analysis/fastq_screen/{sample.sample}_R2_screen.html", sample=samples.itertuples()) +
            expand("analysis/trim_galore/{sample.sample}_R{read}_val_{read}_fastqc.html", sample=samples.itertuples(), read=["1","2"]) +
            expand("analysis/{align_dir}/idxstats/{sample.sample}.idxstats", sample=samples.itertuples(), align_dir=['bowtie2','bowtie2_filt_endogenous']) +
            expand("analysis/{align_dir}/flagstat/{sample.sample}.flagstat", sample=samples.itertuples(), align_dir=['bowtie2','bowtie2_filt_endogenous']) +
            expand("analysis/{align_dir}/CollectAlignmentSummaryMetrics/{sample.sample}.aln_metrics.txt", sample=samples.itertuples(), align_dir=['bowtie2','bowtie2_filt_endogenous']) +
            expand("analysis/{align_dir}/CollectInsertSizeMetrics/{sample.sample}.insert_size_metrics.txt", sample=samples.itertuples(), align_dir=['bowtie2','bowtie2_filt_endogenous']) +
            expand("analysis/{align_dir}/qualimap/{sample.sample}/done", sample=samples.itertuples(), align_dir=['bowtie2','bowtie2_filt_endogenous']) +
            expand("analysis/macs3_narrow/{sample.sample}_peaks.xls", sample=samples.itertuples()) +
            expand("analysis/frags_over_peaks/{peak_type}_frip.csv", peak_type = ['all','group'])
    output:
        mqc_vers="analysis/multiqc/all_mqc_versions.yaml",
        mqc="analysis/multiqc/multiqc_report.html",
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    params:
        config_file=configfile_path,
        workdir="analysis/multiqc/",
        dirs=lambda wildcards,input: " ".join(pd.unique([os.path.dirname(x) for x in input.mqc_dirs])),
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
        --sample-filters {input.sample_filters} \
        --config {input.mqc_config} \
        --outdir {params.workdir} \
        --filename {params.outfile} \
        {params.dirs} {output.mqc_vers}

        """
