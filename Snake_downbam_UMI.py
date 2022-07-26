#configfile: "config.yaml"
# snakemake -j 9 -s Snake_downbam_UMI.py --configfile config/se_downsample.yaml --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn={params.cores} -q home-yeo" -n
# downsample bams to match the # of mapped reads

import pandas as pd
import os
import sys
import glob
import numpy as np

manifest = pd.read_table(config['DOWNSAMPLE_MENIFEST'], index_col = 0, sep = ',')
config['CLIPper_pvalthes'] = None


def get_total_reads(sample_label):
    with open(f"downsample_bam_UMI/{sample_label}.totalcount", 'r') as f:
        nread = int(f.readlines())[0]
    return int(nread.rstrip())

# downsample targets
targets = [int(i) for i in np.array([0.5, 1, 2,  5, 10, 50, 100]) * 10**6]
seeds = []
sample_labels = manifest.Sample.tolist()

module Snake_CLIPper:
    snakefile:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        #"Snake_downbam_UMI",
        "rules/Snake_CLIPper.py"
    config:
        config

module peak_anno:
    snakefile:
        "rules/Snake_peakanno.py"
    config: config


rule all:
    input:
        expand("downsample_bam_UMI_counts/{sample_label}.{nread}.rmDup.count", 
        sample_label = sample_labels,
        nread = targets),
        expand("downsample_bam_UMI/CLIPper/{sample_label}.{nread}.peaks.normed.compressed.annotate.bed", 
        sample_label = ['Dan_singleplex_K562_rep1_RBFOX2', 'katieoligo_RBFOX2_rep1'],
        nread = np.array(targets)[[0,1,2]])
    output:
        "snakeUMI.txt"
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

rule count_read_begin:
    input:
        bam=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["bam"].values[0])
    output:
        "downsample_bam_UMI/{sample_label}.totalcount"
    params:
        run_time=1,
        error_out_file = "error_files/downsample",
        cores = "1",
    shell:
        """
        module load samtools
        samtools view -cF 4 {input.bam} > {output}
        """

rule downsample_bam:
    input:
        bam=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["bam"].values[0]),
        count_csv = "downsample_bam_UMI/{sample_label}.totalcount"
    output:
        subsample_bam="downsample_bam_UMI/{sample_label}.{nread}.bam",
        subsample_bai="downsample_bam_UMI/{sample_label}.{nread}.bam.bai"
    params:
        run_time=1,
        error_out_file = "error_files/downsample",
        cores = "1",
        nread="{nread}"
    shell:
        """
        module load samtools
        total=$(cat {input.count_csv})

        frac=$(awk "BEGIN {{print ({params.nread})/$total}}")
        
        samtools view -h -F 4 -s 5$frac {input.bam} | samtools sort - | samtools view -Sb - > {output.subsample_bam}
        samtools index {output.subsample_bam}
        """

rule umi_dedup:
    input:
        bam="downsample_bam_UMI/{sample_label}.{nread}.bam",
        bai="downsample_bam_UMI/{sample_label}.{nread}.bam.bai"
    output:
        bam_dedup="downsample_bam_UMI/{sample_label}.{nread}.rmDup.bam"
    params:
        error_out_file="error_files/umidedup",
        run_time = 4,
        cores = "4",
        memory = "10000",
        job_name = "sortbam",
        prefix='downsample_bam_UMI/{sample_label}.{nread}'
    shell:
        """
        module load eclip;
        umi_tools dedup \
            --random-seed 1 \
            -I {input.bam} \
            --method unique \
            --output-stats {params.prefix} \
            -S {output.bam_dedup}
        """

rule count_read:
    input:
        bam_dedup="downsample_bam_UMI/{sample_label}.{nread}.rmDup.bam"
    output:
        "downsample_bam_UMI_counts/{sample_label}.{nread}.rmDup.count"
    params:
        error_out_file = "error_files/all",
        run_time = 1,
        cores = "1",
        memory = "20",
        job_name = "all"
    shell:
        """
        module load samtools
        samtools view -cF 4 {input.bam_dedup} > {output}
        """

rule index_bam:
    input:
        bam_dedup="downsample_bam_UMI/{sample_label}.{nread}.rmDup.bam"
    output:
        "downsample_bam_UMI/{sample_label}.{nread}.rmDup.bam.bai"
    params:
        error_out_file = "error_files/all",
        run_time = 1,
        cores = "1",
        memory = "20",
        job_name = "all"
    shell:
        """
        module load samtools
        samtools index {input.bam_dedup} 
        """

use rule clipper from Snake_CLIPper as clipper with:
    input:
        subsample_bam="downsample_bam_UMI/{sample_label}.{nread}.rmDup.bam",
        subsample_bai="downsample_bam_UMI/{sample_label}.{nread}.rmDup.bam.bai"
    output:
        peak="downsample_bam_UMI/CLIPper/{sample_label}.{nread}.peaks.bed"
    params:
        error_out_file = "error_files/all",
        run_time = 16,
        cores = "10",
        memory = "20",
        job_name = "all",
        species=config['SPECIES'],

use rule count_read_num from Snake_CLIPper as count_read2 with:
    input:
        subsample_bam_ip="downsample_bam_UMI/{sample_label}.{nread}.rmDup.bam",
        subsample_bam_in=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["bam_control"].values[0])
        
    output:
        nread_ip="downsample_bam_UMI/output/{sample_label}.{nread}.ip.readnum.txt",
        nread_in="downsample_bam_UMI/output/{sample_label}.{nread}.in.readnum.txt"

use rule norm_peaks from Snake_CLIPper as norm_peak with:
    input:
        subsample_bam_ip="downsample_bam_UMI/{sample_label}.{nread}.rmDup.bam",
        subsample_bam_in=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["bam_control"].values[0]),
        nread_ip="downsample_bam_UMI/output/{sample_label}.{nread}.ip.readnum.txt",
        nread_in="downsample_bam_UMI/output/{sample_label}.{nread}.in.readnum.txt",
        peak="downsample_bam_UMI/CLIPper/{sample_label}.{nread}.peaks.bed"
    output:
        norm_peak="downsample_bam_UMI/CLIPper/{sample_label}.{nread}.peaks.normed.bed"

use rule compress_peak from Snake_CLIPper as compress_peak with:
    input:
        norm_peak="downsample_bam_UMI/CLIPper/{sample_label}.{nread}.peaks.normed.bed"
    output:
        compress_peak="downsample_bam_UMI/CLIPper/{sample_label}.{nread}.peaks.normed.compressed.bed"

use rule annotate from peak_anno as annotate_peak with:
    input:
        peak="downsample_bam_UMI/CLIPper/{sample_label}.{nread}.peaks.normed.compressed.bed"
    output:
        "downsample_bam_UMI/CLIPper/{sample_label}.{nread}.peaks.normed.compressed.annotate.bed"