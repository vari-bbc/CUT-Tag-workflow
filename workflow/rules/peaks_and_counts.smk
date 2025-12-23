rule deeptools_heatmap_peaks:
    input:
        bw=lambda wildcards: expand("analysis/avg_bigwigs/{{scale_method}}/{sample_group}.bw", sample_group=pd.unique(samples_no_controls[samples_no_controls['enriched_factor']==wildcards.enriched_factor]["sample_group"])),
        peaks="analysis/{peak_type}/merged_{enriched_factor}/all.bed",
    output:
        compmat="analysis/deeptools_heatmap_peaks/{scale_method}/{enriched_factor}__{peak_type}_compmat.gz",
        sorted_regions="analysis/deeptools_heatmap_peaks/{scale_method}/{enriched_factor}__{peak_type}_sorted_regions.bed",
        heatmap="analysis/deeptools_heatmap_peaks/{scale_method}/{enriched_factor}__{peak_type}.pdf"
    benchmark:
        "benchmarks/deeptools_heatmap_peaks/{scale_method}/{enriched_factor}__{peak_type}.txt"
    envmodules:
        config['modules']['deeptools']
    params:
        after="2000",
        before="2000",
        binsize=10,
        samp_labels=lambda wildcards, input: " ".join(os.path.basename(x).replace(".bw", "") for x in input.bw),
        temp="analysis/deeptools_heatmap_peaks/{scale_method}/{enriched_factor}__{peak_type}_tmp",
        yaxislabel='"Normalized coverage"',
    threads: 16
    resources:
        mem_gb=100,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        export TMPDIR={params.temp}

        computeMatrix \
        reference-point \
        --referencePoint "center" \
        -p {threads} \
        -b {params.before} \
        -a {params.after} \
        --samplesLabel {params.samp_labels} \
        --binSize {params.binsize} \
        -R {input.peaks} \
        -S {input.bw} \
        -o {output.compmat} 

        echo "END computeMatrix"
        echo "END computeMatrix" 1>&2

        plotHeatmap \
        --heatmapWidth 6 \
        --yAxisLabel {params.yaxislabel} \
        --xAxisLabel "" \
        --refPointLabel "Peak center" \
        --regionsLabel "Peaks" \
        -m {output.compmat} \
        -out {output.heatmap} \
        --outFileSortedRegions {output.sorted_regions}

        echo "END plotHeatmap"
        echo "END plotHeatmap" 1>&2

        """

def get_macs3_bams(wildcards):
    macs3_bams = { 'trt': "analysis/bowtie2_filt_endogenous/{sample}.sorted.bam".format(sample=wildcards.sample) }
    
    control = samples[samples['sample'] == wildcards.sample]['control'].values[0]
    
    if (not pd.isnull(control)):
        macs3_bams['control'] = expand("analysis/bowtie2_filt_endogenous/{sample}.sorted.bam", sample = control.split(','))
    
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
        bed = temp("analysis/frags_over_peaks/{peak_type}/{merge_id}/{sample}.bed"),
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
        bams = lambda wildcards: expand("analysis/bowtie2_filt_{{chromatin_source}}/{sample}.sorted.bam", sample=samples_no_controls[samples_no_controls['enriched_factor'] == wildcards.enriched_factor]['sample']),
    output:
        binned="analysis/csaw_win{width}/csaw_count_{chromatin_source}/{enriched_factor}/binned.rds",
        small_wins="analysis/csaw_win{width}/csaw_count_{chromatin_source}/{enriched_factor}/small_wins.rds",
        filt_small_wins="analysis/csaw_win{width}/csaw_count_{chromatin_source}/{enriched_factor}/filt_small_wins.rds",
        global_filt="analysis/csaw_win{width}/csaw_count_{chromatin_source}/{enriched_factor}/global_filt.rds",
        bkgrd_scale="analysis/csaw_win{width}/csaw_count_{chromatin_source}/{enriched_factor}/bkgrd_norm_factors.tsv",
        hiAbund_scale="analysis/csaw_win{width}/csaw_count_{chromatin_source}/{enriched_factor}/hiAbund_norm_factors.tsv",
        bkgrd_scale_link="analysis/bigwig_norm_factors/{enriched_factor}_csaw_win{width}_bkgd_{chromatin_source}.tsv",
        hiAbund_scale_link="analysis/bigwig_norm_factors/{enriched_factor}_csaw_win{width}_hiAbund_{chromatin_source}.tsv",
        figs_dir=directory("analysis/csaw_win{width}/csaw_count_{chromatin_source}/{enriched_factor}/figures")
    benchmark:
        "benchmarks/csaw_win{width}/csaw_count_{chromatin_source}/{enriched_factor}.txt"
    params:
        window_width="{width}",
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
        filt_small_wins="analysis/csaw_win{width}/csaw_count_endogenous/{enriched_factor}/filt_small_wins.rds",
        bin_counts="analysis/csaw_win{width}/csaw_count_endogenous/{enriched_factor}/binned.rds",
        filt_small_wins_spikein="analysis/csaw_win{width}/csaw_count_spikein/{enriched_factor}/filt_small_wins.rds",
        bin_counts_spikein="analysis/csaw_win{width}/csaw_count_spikein/{enriched_factor}/binned.rds",
        merged_peaks="analysis/{peak_type}/merged_{enriched_factor}/all.bed",
        contrasts="bin/contrasts.tsv"
    output:
        edgeR_rds="analysis/csaw_win{width}/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/edgeR_objs.rds",
        results_rds="analysis/csaw_win{width}/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/results.rds",
        results_se_rds="analysis/csaw_win{width}/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/results_se.rds",
        peaks_res_rds="analysis/csaw_win{width}/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/peaks_res.rds",
        combined_rds="analysis/csaw_win{width}/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/combined.rds",
        peaks_anno_rds="analysis/csaw_win{width}/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/peaks_anno.rds",
        figs_dir=directory("analysis/csaw_win{width}/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/figures")
    benchmark:
        "benchmarks/csaw_win{width}/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}.txt"
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
        edgeR_rds="analysis/csaw_win{width}/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/edgeR_objs.rds",
        peaks_res_rds="analysis/csaw_win{width}/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/peaks_res.rds",
        combined_rds="analysis/csaw_win{width}/csaw_diff/{enriched_factor}__{peak_type}__{norm_type}/combined.rds",
    output:
        rmd="analysis/csaw_win{width}/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/csaw_summary.Rmd",
        edgeR_rds="analysis/csaw_win{width}/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/edgeR_objs.rds",
        peaks_res_rds="analysis/csaw_win{width}/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/peaks_res.rds",
        combined_rds="analysis/csaw_win{width}/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/combined.rds",
        out_res="analysis/csaw_win{width}/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/DB_results.xlsx",
        html_report="analysis/csaw_win{width}/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/csaw_summary.html",
        figdir=directory("analysis/csaw_win{width}/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}/individual_figures")
    benchmark:
        "benchmarks/csaw_win{width}/csaw_summary/{enriched_factor}__{peak_type}__{norm_type}.txt"
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

