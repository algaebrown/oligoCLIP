# simulate at different sequencing depth, how many UMI/peaks will be captured (single end)
import pandas as pd
import glob
import numpy as np

print(config)
manifest = pd.read_table(config['DOWNSAMPLE_MENIFEST'], index_col = 0, sep = ',')

#snakemake -j 12 -s Snake_downbam_UMI_PE.py --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn={params.cores} -q home-yeo" --configfile config/pe_downsample.yaml

def get_total_reads(sample_label):
    with open(f"downsample_bam_UMI/{sample_label}.totalcount", 'r') as f:
        nread = int(f.readlines())[0]
    return int(nread.rstrip())

# downsample targets
targets = [int(i) for i in np.array([1, 2, 5, 10, 50, 100]) * 10**6]
seeds = []
sample_labels = manifest.Sample.tolist()


module Snake_downbam:
    snakefile:
        "Snake_downbam_UMI.py"
    config: config
        
module Snake_CLIPper:
    snakefile:
        "Snake_CLIPper.py"
    config: config

module peak_anno:
    snakefile:
        "Snake_peakanno.py"
    config: config

combinations = []
for index, row in manifest.iterrows():
    combinations.append(row['Sample']+'-'+row['barcode'])

rule all:
    input:
        expand("downsample_bam_UMI_counts/{combination}.{nread}.rmDup.count", combination=combinations,nread = targets)+
        expand("downsample_bam_UMI/CLIPper/{sample_label}.{nread}.peaks.normed.compressed.annotate.bed", 
        sample_label = ['676_01_RBFOX2'],
        nread = np.array(targets)[[0,1,2]]
        )
    output:
        "snakeUMI_PE.txt"
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

use rule count_read_begin from Snake_downbam as pe_countread with:
    input:
        bam=lambda wildcards: glob.glob(manifest.loc[
                (manifest.Sample == wildcards.sample_label)&
                (manifest.barcode == wildcards.barcode)
            ]["bam"].values[0])
    output:
        "downsample_bam_UMI/{sample_label}-{barcode}.totalcount"

use rule downsample_bam from Snake_downbam as PE_downsample with:
    input:
        bam=lambda wildcards: glob.glob(manifest.loc[
                    (manifest.Sample == wildcards.sample_label)&
                    (manifest.barcode == wildcards.barcode)
                ]["bam"].values[0]),
        count_csv = "downsample_bam_UMI/{sample_label}-{barcode}.totalcount"
    output:
        subsample_bam="downsample_bam_UMI/{sample_label}-{barcode}.{nread}.bam",
        subsample_bai="downsample_bam_UMI/{sample_label}-{barcode}.{nread}.bam.bai"
    params:
        nread = lambda wildcards: int(wildcards.nread)/(manifest.loc[manifest.Sample == wildcards.sample_label].shape[0]),
        run_time=1,
        error_out_file = "error_files/downsample",
        cores = "1",

rule PE_dedup:
    input:
        bam="downsample_bam_UMI/{sample_label}-{barcode}.{nread}.bam",
        bai="downsample_bam_UMI/{sample_label}-{barcode}.{nread}.bam.bai"
    output:
        bam_dedup="downsample_bam_UMI/{sample_label}-{barcode}.{nread}.rmDup.bam"
    params:
        error_out_file="error_files/umidedup",
        run_time = 4,
        cores = "4",
        memory = "10000",
        job_name = "sortbam",
        prefix='downsample_bam_UMI/{sample_label}.{barcode}.{nread}'
    shell:
        """
        module load eclip;
        umi_tools dedup \
            --random-seed 1 \
            -I {input.bam} \
            --method unique \
            --output-stats {params.prefix} \
            --paired \
            -S {output.bam_dedup}
        """


rule merge_bams:
    input:
        lambda wildcards: expand("downsample_bam_UMI/{sample_label}-{barcode}.{nread}.rmDup.bam", 
            sample_label = [wildcards.sample_label],
            nread = [wildcards.nread],
            barcode = manifest.loc[manifest.Sample == wildcards.sample_label, 'barcode'].tolist()
        )
    output:
        "downsample_bam_UMI/{sample_label}.{nread}.rmDup.merged.bam"
    params:
        error_out_file="error_files/umidedup",
        run_time = 2,
        cores = "4",
        memory = "10000",
        job_name = "sortbam",
    shell:
        """
        module load samtools;
        samtools merge {output} {input}
        """

rule sort_and_index_extract_r2:
    input:
        "downsample_bam_UMI/{sample_label}.{nread}.rmDup.merged.bam"
    output:
        bam = "downsample_bam_UMI/{sample_label}.{nread}.rmDup.merged.r2.bam",
        bai="downsample_bam_UMI/{sample_label}.{nread}.rmDup.merged.r2.bam.bai"
    params:
        error_out_file="error_files/umidedup",
        run_time = 1,
        cores = "4",
        memory = "10000",
        job_name = "sortmerge",
    shell:
        """
        module load samtools;
        samtools view -h -f 0x0080 {input} | samtools sort - | samtools view -Sb - > {output.bam}
        samtools index {output.bam}
        """
rule count_read:
    input:
        bam_dedup="downsample_bam_UMI/{sample_label}-{barcode}.{nread}.rmDup.bam",
    output:
        "downsample_bam_UMI_counts/{sample_label}-{barcode}.{nread}.rmDup.count"
    params:
        error_out_file = "error_files/all",
        run_time = 1,
        cores = "1",
        memory = "20",
        job_name = "all"
    shell:
        """
        module load samtools
        samtools view -cF 4 -f 0x80 {input.bam_dedup} > {output}
        """

use rule clipper from Snake_CLIPper as PE_clipper with:
    input:
        subsample_bam="downsample_bam_UMI/{sample_label}.{nread}.rmDup.merged.r2.bam",
        subsample_bai="downsample_bam_UMI/{sample_label}.{nread}.rmDup.merged.r2.bam.bai"
    output:
        peak="downsample_bam_UMI/CLIPper/{sample_label}.{nread}.peaks.bed"
    params:
        error_out_file = "error_files/all",
        run_time = 16,
        cores = "10",
        memory = "20",
        job_name = "all",
        species=config['SPECIES']

use rule count_read_num from Snake_CLIPper as count_read2 with:
    input:
        subsample_bam_ip="downsample_bam_UMI/{sample_label}.{nread}.rmDup.merged.r2.bam",
        subsample_bam_in=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["bam_control"].values[0])
        
    output:
        nread_ip="downsample_bam_UMI/output/{sample_label}.{nread}.ip.readnum.txt",
        nread_in="downsample_bam_UMI/output/{sample_label}.{nread}.in.readnum.txt"

use rule norm_peaks from Snake_CLIPper as norm_peak with:
    input:
        subsample_bam_ip="downsample_bam_UMI/{sample_label}.{nread}.rmDup.merged.r2.bam",
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