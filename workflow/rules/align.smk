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
        sorted_bam=temp("analysis/bowtie2/{sample}.sorted.bam") if config['cleanup_big_files']['unfilt_coordsorted_bams'] else "analysis/bowtie2/{sample}.sorted.bam",
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
        samtools sort -m 6G -O "bam" -@ {threads} -o {output.sorted_bam} -
        
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
            grep -P "^$chr\\t" {output.keep_regions_temp} >> {output.keep_regions} 
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
        sorted_bam=temp("analysis/{align_dirname}_filt_endogenous/{bam_name}.sorted.bam") if config['cleanup_big_files']['filt_coordsorted_bams'] else "analysis/{align_dirname}_filt/{bam_name}.sorted.bam",
        sorted_bai="analysis/{align_dirname}_filt_endogenous/{bam_name}.sorted.bam.bai",
        bam=temp("analysis/{align_dirname}_filt_endogenous/{bam_name}.bam") if config['cleanup_big_files']['filt_namecollated_bams'] else "analysis/{align_dirname}_filt/{bam_name}.bam",
    params:
        view_mapq="" if config['samtools_mapq'] == "" else "-q {mapq}".format(mapq=config['samtools_mapq']),
        view_keep="" if config['samtools_keep_flags'] == "" else "-f {flag}".format(flag=config['samtools_keep_flags']),
        view_omit="" if config['samtools_omit_flags'] == "" else "-F {flag}".format(flag=config['samtools_omit_flags']),
    benchmark:
        "benchmarks/{align_dirname}_filt_endogenous/{bam_name}.txt"
    envmodules:
        config['modules']['samtools']
    threads: 8
    resources:
        mem_gb = 96,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        samtools view -@ {threads} -u --region-file {input.keep_regions} {params.view_mapq} {params.view_keep} {params.view_omit} {input.bam} | \
        samtools collate -@ {threads} -O -u - | \
        samtools fixmate -@ {threads} -u - - | \
        samtools view -@ {threads} -f 2 -b -u - | tee {output.bam} | \
        samtools sort -m 6G -@ {threads} -o {output.sorted_bam} -

        samtools index -@ {threads} {output.sorted_bam}

        """


rule get_spikein_chroms_bed:
    input:
        ref_fai=config["ref"]["fai"],
        spikein_chroms = config['spikein_chroms_file']
    output:
        keep_fai=temp("analysis/misc/spikein_chroms.fai"),
        keep_bed="analysis/misc/spikein_chroms.bed"
    benchmark:
        "benchmarks/get_spikein_chroms_bed/bench.txt"
    params:
    threads: 1
    resources:
        mem_gb=20,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
    shell:
        """
        for chr in $(cat {input.spikein_chroms})
        do
            grep -P "^$chr\\t" {input.ref_fai} >> {output.keep_fai}
        done

        perl -lane 'print qq|$F[0]\\t0\\t$F[1]|' {output.keep_fai} > {output.keep_bed}
        """

rule filter_spikein_bams:
    """
    Filter BAMs for spikeins. After filters, fixmate will be run to update SAM flags for pairs where one mate is removed and then filtered for properly paired alignments
    """
    input:
        bam = "analysis/{align_dirname}/{bam_name}.sorted.bam",
        spikein_chroms = "analysis/misc/spikein_chroms.bed",
    output:
        sorted_bam="analysis/{align_dirname}_filt_spikein/{bam_name}.sorted.bam",
        sorted_bai="analysis/{align_dirname}_filt_spikein/{bam_name}.sorted.bam.bai",
        bam="analysis/{align_dirname}_filt_spikein/{bam_name}.bam",
    params:
        view_mapq="" if config['samtools_mapq'] == "" else "-q {mapq}".format(mapq=config['samtools_mapq']),
        view_keep="" if config['samtools_keep_flags'] == "" else "-f {flag}".format(flag=config['samtools_keep_flags']),
        view_omit="" if config['samtools_omit_flags'] == "" else "-F {flag}".format(flag=config['samtools_omit_flags']),
    benchmark:
        "benchmarks/{align_dirname}_filt_spikein/{bam_name}.txt"
    envmodules:
        config['modules']['samtools']
    threads: 8
    resources:
        mem_gb = 96,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        samtools view -@ {threads} -u --region-file {input.spikein_chroms} {params.view_mapq} {params.view_keep} {params.view_omit} {input.bam} | \
        samtools collate -@ {threads} -O -u - | \
        samtools fixmate -@ {threads} -u - - | \
        samtools view -@ {threads} -f 2 -b -u - | tee {output.bam} | \
        samtools sort -m 6G -@ {threads} -o {output.sorted_bam} -

        samtools index -@ {threads} {output.sorted_bam}

        """