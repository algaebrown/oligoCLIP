
import pandas as pd
#snakemake -j 8 -s Snake_peakanno_chi.py --configfile config/annotate_chi/multiplex.yaml --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn={params.cores} -q home-yeo" -n --use-conda

try:
    WORKDIR=config['WORKDIR']
    workdir: WORKDIR

    print('WORKDIR', os.getcwd())
    sample_labels,= glob_wildcards("input/{sample_label}.peak.chisq.normed.bed")
    
    
    GTF = config['GTF']
    SCRIPT_PATH=config['SCRIPT_PATH']
    ANNOTATOR_SPECIES = config['ANNOTATOR_SPECIES']
    GENOMEFA = config['GENOMEFA']
    CHI_LOG10PVAL_THRES = config['CHI_LOG10PVAL_THRES']
    config['CLIPper_pvalthes'] = None # to use module

    all_rbfox = [s for s in sample_labels if 'RBFOX' in s]
except Exception as e:
    print('config:', e)


module motif:
    snakefile:
        "rules/snake_scoreRBNS_SELEX.py"
    config: config

module peak_anno:
    snakefile:
        "rules/Snake_peakanno.py"
    config: config



rule all:
    input:
        expand("output/{sample_label}.peak.chisq.normed.annotate.bed", sample_label = sample_labels),
        # expand("output/{sample_label}peak.chisq.normed.  gc", sample_label = sample_labels)+
        # expand("output/{sample_label}peak.chisq.normed.  filtered.bed", sample_label = sample_labels)+
        expand("output/motif/{sample_label}.peak.chisq.normed.filtered.annotate.svg", sample_label = sample_labels),
        expand("output/{sample_label}.peak.chisq.normed.annotate.motifscore.csv", sample_label = all_rbfox),
        expand("summary/chisq/{sample_label}.peaks.summary", sample_label = sample_labels)
    output:
        "snakepeak_annochi.txt"
    params:
        error_out_file = "error_files/all",
        run_time = 1,
        cores = "1",
        memory = "20",
        job_name = "all"
    shell:
        """
        echo $(date)  > {output};
        echo created by HLH and the Yeo lab >> {output}
        """


use rule annotate from peak_anno as annotate with:
    input:
        peak="input/{sample_label}.peak.chisq.normed.bed"
    output:
        "output/{sample_label}.peak.chisq.normed.annotate.bed"


use rule calc_partition_nuc from peak_anno as calc_partition_nuc with:
    input:
        partition = "input/{sample_label}.peak.chisq.normed.bed",
        genome = GENOMEFA
    output:
        nuc = "output/{sample_label}.peak.chisq.normed.nuc",
        gc = "output/{sample_label}.peak.chisq.normed.gc",


rule filter_peak_fc_pval:
    input:
        "input/{sample_label}.peak.chisq.normed.bed"
    output:
        "output/{sample_label}.peak.chisq.normed.filtered.bed"
    params:
        run_time=1,
        error_out_file = "error_files/annotate",
        cores = "1",
        thres = CHI_LOG10PVAL_THRES, # log10 pval
    shell: 
        """
        awk '{{ if ($5 > {params.thres}) print }}' {input}  > {output}
        """


use rule motif_analysis from peak_anno as motif_analysis with:
    input:
        peak="output/{sample_label}.peak.chisq.normed.filtered.bed"
    output:
        pickle="output/motif/{sample_label}.peak.chisq.normed.filtered.annotate.pickle",
        svg="output/motif/{sample_label}.peak.chisq.normed.filtered.annotate.svg",
    params:
        run_time=12,
        error_out_file = "error_files/annotate",
        cores = "1",
        fa = config['GENOMEFA'],
        sps = config['ANNOTATOR_SPECIES'],
        homer="output/CLIPper/{sample_label}.peaks.chisq.normed.filtered.annotate.homer",

rule peak_summary_chi:
    input:
        peak="output/{sample_label}.peak.chisq.normed.annotate.bed"
    output:
        summary_csv = "summary/chisq/{sample_label}.peaks.summary"
    params:
        run_time=1,
        error_out_file = "error_files/annotate",
        cores = "1",
    conda:
        "rules/envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/summarize_peak.py --bed {input.peak}  \
                                        --is_chi \
                                        --outfile {output.summary_csv}
        """
    
use rule * from motif as motif_*