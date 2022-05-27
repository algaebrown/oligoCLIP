#DOWNSAMPLE to match the # of UMI
# snakemake -j 70 -s Snake_downbam.py --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn=1 -q home-yeo -e {params.error_out_file} -o /dev/null"
#snakemake -j 40 -s Snake_downbam.py --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn=1 -q home-yeo -e {params.error_out_file} -o /dev/null" downsample_bams/{ENCODE_Dan_multiplex1_K562_rep3_IGF2BP2_CLIP,Dan_multiplex1_K562_rep3_IGF2BP2_down_INPUT,ENCODE_Dan_multiplex1_K562_rep5_IGF2BP2_CLIP,Dan_multiplex1_K562_rep5_IGF2BP2_down_INPUT,ENCODE_Dan_multiplex1_K562_rep3_RBFOX2_CLIP,Dan_multiplex1_K562_rep3_RBFOX2_down_INPUT,ENCODE_Dan_multiplex1_K562_rep5_RBFOX2_CLIP,Dan_multiplex1_K562_rep5_RBFOX2_down_INPUT,ENCODE_Dan_multiplex1_K562_rep3_PUM2_CLIP,Dan_multiplex1_K562_rep3_PUM2_down_INPUT,ENCODE_Dan_multiplex1_K562_rep5_PUM2_CLIP,Dan_multiplex1_K562_rep5_PUM2_down_INPUT,Dan_multiplex1_K562_rep3_FAM120A_down_CLIP,Dan_multiplex1_K562_rep3_FAM120A_down_INPUT,ENCODE_Dan_multiplex1_K562_rep5_FAM120A_CLIP,Dan_multiplex1_K562_rep5_FAM120A_down_INPUT,ENCODE_Dan_multiplex1_K562_rep3_DDX3_CLIP,Dan_multiplex1_K562_rep3_DDX3_down_INPUT,ENCODE_Dan_multiplex1_K562_rep5_DDX3_CLIP,Dan_multiplex1_K562_rep5_DDX3_down_INPUT,ENCODE_Dan_multiplex1_K562_rep3_ZC3H11A_CLIP,Dan_multiplex1_K562_rep3_ZC3H11A_down_INPUT,ENCODE_Dan_multiplex1_K562_rep5_ZC3H11A_CLIP,Dan_multiplex1_K562_rep5_ZC3H11A_down_INPUT,ENCODE_Dan_multiplex1_K562_rep3_EIF3G_CLIP,Dan_multiplex1_K562_rep3_EIF3G_down_INPUT,ENCODE_Dan_multiplex1_K562_rep5_EIF3G_CLIP,Dan_multiplex1_K562_rep5_EIF3G_down_INPUT,ENCODE_Dan_multiplex1_K562_rep3_PRPF8_CLIP,Dan_multiplex1_K562_rep3_PRPF8_down_INPUT,ENCODE_Dan_multiplex1_K562_rep5_PRPF8_CLIP,Dan_multiplex1_K562_rep5_PRPF8_down_INPUT,ENCODE_Dan_multiplex1_K562_rep3_LIN28B_CLIP,Dan_multiplex1_K562_rep3_LIN28B_down_INPUT,ENCODE_Dan_multiplex1_K562_rep5_LIN28B_CLIP,Dan_multiplex1_K562_rep5_LIN28B_down_INPUT,ENCODE_Dan_multiplex1_K562_rep3_SF3B4_CLIP,Dan_multiplex1_K562_rep3_SF3B4_down_INPUT,ENCODE_Dan_multiplex1_K562_rep5_SF3B4_CLIP,Dan_multiplex1_K562_rep5_SF3B4_down_INPUT}.bam -F
# downsample bams to match the # of mapped reads

import pandas as pd
import os
import glob


MANIFEST='/home/hsher/projects/oligo_results/downsample_menifest.csv' 
# columns: bam, target_nread, mapped_nread, Sample
# bam: filename to the bam file (deduplicated)
# target_nread: number of read to downsample to
# mapped_nread: number of read the bam file contains
# Sample: sample name, needs to be unique

if not os.path.exists(MANIFEST): make_meta(MANIFEST)
manifest = pd.read_table(MANIFEST, index_col = 0, sep = ',')
print(manifest)

manifest['frac']=manifest['target_nread']/manifest['mapped_nread']
sample_labels = manifest.Sample.tolist()

rule all:
    input:
        expand("downsample_bams/{sample_label}.bam.bai", sample_label = sample_labels)
        
    output:
        "snakeCLIP.txt"
    params:
        error_out_file = "error_files/all",
        run_time = "0:10:00",
        cores = "1",
        memory = "20",
        job_name = "all"
    shell:
        """
        echo $(date)  > {output};
        echo created by HLH and the Yeo lab >> {output}
        """

rule downsample_bam:
    input:
        bam=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["bam"].values[0]),
    output:
        subsample_bam="downsample_bams/{sample_label}.bam",
        subsample_bai="downsample_bams/{sample_label}.bam.bai"
    params:
        run_time="1:00:00",
        error_out_file = "error_files/downsample",
        frac=lambda wildcards: manifest.loc[manifest.Sample == wildcards.sample_label]["frac"].values[0], # least reads in all 3 IPs
        cores = "1",
    shell:
        """
        module load samtools
        samtools view -h -F 4 -s 5{params.frac} {input.bam} | samtools view -Sb - > {output.subsample_bam}
        samtools index {output.subsample_bam}
        """