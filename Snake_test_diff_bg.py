# call peaks and normalize peaks lah
#configfile: "config.yaml"

#snakemake -j 12 -s Snake_test_diff_bg.py --configfile config/test_INPUT/SLBP.yaml --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn={params.cores} -q home-yeo" 
# downsample IP.bam
# run CLIPper again
# normalize to background again
from gettext import find
from tabnanny import check
from turtle import back
import pandas as pd
import os
import sys
import glob

workdir: config['WORKDIR']
SCRIPT_PATH = config['SCRIPT_PATH']
print(config)

clipper_peaks = {'katieoligo_SLBP_rep1': '/home/hsher/scratch/katie_drosphila/output/CLIPper/katieoligo_SLBP_rep1.peaks.bed',
                'katieoligo_SLBP_rep2': '/home/hsher/scratch/katie_drosphila/output/CLIPper/katieoligo_SLBP_rep2.peaks.bed',
                'Dan_singleplex_K562_rep1_SLBP': '/home/hsher/scratch/downsample_peak_from_eclipse/output/CLIPper/Dan_singleplex_K562_rep1_SLBP.peaks.bed',
                'Dan_singleplex_K562_rep2_SLBP': '/home/hsher/scratch/downsample_peak_from_eclipse/output/CLIPper/Dan_singleplex_K562_rep2_SLBP.peaks.bed'}

ip_bams = {'katieoligo_SLBP_rep1': '/home/hsher/scratch/katie_drosphila/bams/genome/SLBP_CLIP1.genome-mappedSoSo.rmDupSo.Aligned.out.bam',
                'katieoligo_SLBP_rep2': '/home/hsher/scratch/katie_drosphila/bams/genome/SLBP_CLIP2.genome-mappedSoSo.rmDupSo.Aligned.out.bam',
                'Dan_singleplex_K562_rep1_SLBP': '/home/hsher/scratch/ABC_reprocess/K562_SLBP_rep1/bams/genome/SLBP.genome-mappedSoSo.rmDupSo.Aligned.out.bam',
                'Dan_singleplex_K562_rep2_SLBP': '/home/hsher/scratch/ABC_reprocess/K562_SLBP_rep2/bams/genome/SLBP.genome-mappedSoSo.rmDupSo.Aligned.out.bam'}


k562_backgrounds = {'SLBP_SMINPUT':'/projects/ps-yeolab3/encode/analysis/encode_GRCh38_v1/262_INPUT_ATTCAGAA-TATAGCCT_L006_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.bam',
                    'PRPF8_SMINPUT': '/home/hsher/scratch/encode_recall/encode3.278_01_PRPF8_INPUT.NIL.r1.fq.genome-mappedSo.rmDupSo.r2.bam',
                    'K562_multiplex_rep4': '/home/hsher/scratch/k562_rep4.sorted.bam',
                    }
hek_backgrounds = {'HEK293T_RNAseq': '/home/hsher/seqdata/20210702_293t_ribodep_rnaseq/SRR8181090Aligned.sortedByCoord.out.bam',
                    'SLBP_SMINPUT': '/home/hsher/scratch/SLBP/SLBP_nature2016/results/Nature2016_SLBP.SLBP_13_INPUT.NIL.r1.fq.genome-mappedSo.rmDupSo.r2.bam',
                    'RPS2_INPUT': '/home/hsher/scratch/katie_drosphila/bams/genome/RPS2_INPUT1.genome-mappedSoSo.rmDupSo.Aligned.out.bam'}

background_sets = {'katieoligo_SLBP_rep1': hek_backgrounds,
                    'katieoligo_SLBP_rep2': hek_backgrounds,
                'Dan_singleplex_K562_rep1_SLBP': k562_backgrounds,
                'Dan_singleplex_K562_rep2_SLBP': k562_backgrounds
}



def check_path(d):
    for i in d.values():
        if not os.path.isfile(i):
            print(i)
        
check_path(clipper_peaks)
check_path(ip_bams)
check_path(k562_backgrounds)
check_path(hek_backgrounds)

module peak_anno:
    snakefile:
        "rules/Snake_peakanno.py"
    config: config

module peak_call:
    snakefile:
        "rules/Snake_CLIPper.py"
    config: config

def find_all_output(wildcards):
    all_prefix = []
    for ip in clipper_peaks.keys():
        for in_ in background_sets[ip].keys():
            #combination = f"output/{ip}/{in_}.peaks.normed.compressed.annotate.bed"
            combination = f"summary/normalized/{ip}/{in_}.peaks.summary"
            all_prefix.append(combination)
    return all_prefix

print(find_all_output(config))

rule all:
    input:
        lambda wildcards: find_all_output(wildcards)
    output:
        "snakeCLIPtestinput.txt"
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

rule count_read_num:
    input:
        bam_in=lambda wildcards: glob.glob(background_sets[wildcards.ip_name][wildcards.in_name])
    output:
        nread_in="output/{ip_name}/{in_name}.in.readnum.txt",
    params:
        run_time=2,
        error_out_file = "error_files/countread",
        cores = "1",
    shell:
        """
        module load samtools;
        samtools view -cF 4 {input.bam_in} > {output.nread_in}
        """

use rule norm_peaks from peak_call as ssnorm_peak with:
    input:
        subsample_bam_ip=lambda wildcards: ip_bams[wildcards.ip_name],
        subsample_bam_in=lambda wildcards: background_sets[wildcards.ip_name][wildcards.in_name],
        nread_ip=lambda wildcards: clipper_peaks[wildcards.ip_name].replace('.peaks.bed', '.ip.readnum.txt').replace('CLIPper/', ''),
        nread_in="output/{ip_name}/{in_name}.in.readnum.txt",
        peak=lambda wildcards: clipper_peaks[wildcards.ip_name],
    output:
        norm_peak="output/{ip_name}/{in_name}.peaks.normed.bed"

use rule compress_peak from peak_call as compress_peak with:
    input:
        norm_peak="output/{ip_name}/{in_name}.peaks.normed.bed"
    output:
        compress_peak="output/{ip_name}/{in_name}.peaks.normed.compressed.bed"

use rule annotate from peak_anno as annotate_peak with:
    input:
        peak="output/{ip_name}/{in_name}.peaks.normed.compressed.bed"
    output:
        "output/{ip_name}/{in_name}.peaks.normed.compressed.annotate.bed"
use rule peak_summary_normalized from peak_anno with:
    input:
        peak="output/{ip_name}/{in_name}.peaks.normed.compressed.annotate.bed"
    output:
        summary_csv = "summary/normalized/{ip_name}/{in_name}.peaks.summary"
        