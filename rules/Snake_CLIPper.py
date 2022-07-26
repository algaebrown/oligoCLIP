# some default


import glob
import pandas as pd


try:
    manifest = pd.read_table(config['MANIFEST'], index_col = False, sep = ',')
except Exception as e:
    print(e)
    pass

rule clipper:
    input:
        bam=lambda wildcards: glob.glob(manifest.loc[manifest.uid == wildcards.sample_label]["bam_0"].values[0]),
        bai=lambda wildcards: glob.glob(manifest.loc[manifest.uid == wildcards.sample_label]["bam_0"].values[0]+'.bai')
    output:
        peak="output/CLIPper/{sample_label}.peaks.bed"
        
    params:
        run_time=16,
        species=config['SPECIES'],
        error_out_file = "error_files/clipper",
        cores = "10",
    shell:
        """
        module load clipper/42502ec
        clipper -b {input.bam} -s {params.species} -o {output.peak} --processors=16
        """
rule count_read_num:
    input:
        subsample_bam_ip=lambda wildcards: glob.glob(manifest.loc[manifest.uid == wildcards.sample_label]["bam_0"].values[0]),
        subsample_bam_in=lambda wildcards: glob.glob(manifest.loc[manifest.uid == wildcards.sample_label]["bam_control_0"].values[0])
        
    output:
        nread_ip="output/{sample_label}.ip.readnum.txt",
        nread_in="output/{sample_label}.in.readnum.txt"
        
    params:
        run_time=2,
        error_out_file = "error_files/countread",
        cores = "1",

    shell:
        """
        module load samtools;
        samtools view -cF 4 {input.subsample_bam_ip} > {output.nread_ip};
        samtools view -cF 4 {input.subsample_bam_in} > {output.nread_in}
        """
rule norm_peaks:
    input:
        subsample_bam_ip=lambda wildcards: glob.glob(manifest.loc[manifest.uid == wildcards.sample_label]["bam_0"].values[0]),
        subsample_bam_in=lambda wildcards: glob.glob(manifest.loc[manifest.uid == wildcards.sample_label]["bam_control_0"].values[0]),
        nread_ip="output/{sample_label}.ip.readnum.txt",
        nread_in="output/{sample_label}.in.readnum.txt",
        peak="output/CLIPper/{sample_label}.peaks.bed"
        #peak = lambda wildcards: glob.glob(manifest.loc[manifest.uid == wildcards.sample_label]["clipper"].values[0]),
    output:
        norm_peak="output/{sample_label}.peaks.normed.bed"
    params:
        run_time=6,
        error_out_file = "error_files/norm_peaks",
        cores = "4",
        script_path=os.path.join(config['SCRIPT_PATH'], 'overlap_peakfi_with_bam.pl')
    shell:
        """
        module load eclip/0.7.0;
        perl {params.script_path} \
            {input.subsample_bam_ip} \
            {input.subsample_bam_in} \
            {input.peak} \
            {input.nread_ip} \
            {input.nread_in} \
            {output.norm_peak}
        """
rule compress_peak:
    input:
        norm_peak="output/{sample_label}.peaks.normed.bed"
    output:
        compress_peak="output/{sample_label}.peaks.normed.compressed.bed"
    params:
        run_time=14,
        error_out_file = "error_files/compress_peak",
        cores = "1",
    shell:
        """
        module load eclip/0.7.0;
        perl $ECLIP_HOME/bin/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.pl \
        {input.norm_peak} \
        {output.compress_peak}
        """