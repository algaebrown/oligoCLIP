import glob
import pandas as pd
try:
    manifest = pd.read_table(config['MANIFEST'], index_col = False, sep = ',')
except Exception as e:
    print(e)
    pass
def get_all_other_bams(sample_label, lib):
    return manifest.loc[(manifest['lib'] == lib) & (manifest['uid'] != sample_label), 'bam_0'].tolist()

print([len(get_all_other_bams(s, l)) for s, l in zip(manifest['uid'], manifest['lib'])])

module peak_anno:
    snakefile:
        "Snake_peakanno.py"
    config: config

module peak_call:
    snakefile:
        "Snake_CLIPper.py"
    config: config

rule concat_other_bams:
    input:
        other_bams = lambda wildcards: get_all_other_bams(wildcards.sample_label, wildcards.lib)
    output:
        combined_bam = temp("chisq/bg/{lib}/{sample_label}.bam"),
        combined_sorted_bam = "chisq/bg/{lib}/{sample_label}.sorted.bam",
        combined_sorted_index = "chisq/bg/{lib}/{sample_label}.sorted.bam.bai"
    params:
        run_time=2,
        error_out_file = "error_files/annotate",
        cores = "1",
    shell:
        """
        module load samtools
        samtools merge {output.combined_bam} {input.other_bams}
        samtools sort {output.combined_bam} | samtools view -Sb - > {output.combined_sorted_bam}
        samtools index {output.combined_sorted_bam}
        """

use rule count_read_num from peak_call as count_bam_chi with:
    input:
        subsample_bam_ip=lambda wildcards: glob.glob(manifest.loc[manifest.uid == wildcards.sample_label]["bam_0"].values[0]),
        subsample_bam_in="chisq/bg/{lib}/{sample_label}.sorted.bam"
    output:
        nread_ip="chisq/counts/{lib}/{sample_label}.ip.readnum.txt",
        nread_in="chisq/counts/{lib}/{sample_label}.in.readnum.txt"

use rule norm_peaks from peak_call as norm_chisq_bg with:
    input:
        subsample_bam_ip=lambda wildcards: glob.glob(manifest.loc[manifest.uid == wildcards.sample_label]["bam_0"].values[0]),
        subsample_bam_in="chisq/bg/{lib}/{sample_label}.sorted.bam",
        nread_ip="chisq/counts/{lib}/{sample_label}.ip.readnum.txt",
        nread_in="chisq/counts/{lib}/{sample_label}.in.readnum.txt",
        peak="output/CLIPper/{sample_label}.peaks.bed"
    output:
        norm_peak = "chisq/{lib}/{sample_label}.peaks.normed.bed"

use rule compress_peak from peak_call as compress_chi with:
    input:
        norm_peak = "chisq/{lib}/{sample_label}.peaks.normed.bed"
    output:
        compress_peak = "chisq/{lib}/{sample_label}.peaks.normed.compressed.bed"

use rule annotate from peak_anno as anno_chi with:
    input:
        peak="chisq/{lib}/{sample_label}.peaks.normed.compressed.bed"
    output:
        "chisq/{lib}/{sample_label}.peaks.normed.compressed.annotate.bed"

use rule filter_peak_fc_pval from peak_anno as filter_chi with:
    input:
        "chisq/{lib}/{sample_label}.peaks.normed.compressed.bed"
    output:
        "chisq/{lib}/{sample_label}.peaks.normed.compressed.filtered.bed"

use rule motif_analysis from peak_anno as chi_motif with:
    input:
        peak="chisq/{lib}/{sample_label}.peaks.normed.compressed.filtered.bed"
    output:
        pickle="chisq/motif/{lib}/{sample_label}.peaks.normed.compressed.filtered.annotate.pickle",
        svg="chisq/motif/{lib}/{sample_label}.peaks.normed.compressed.filtered.annotate.svg",
    params:
        run_time=3,
        error_out_file = "error_files/annotate",
        cores = "1",
        fa = config['GENOMEFA'],
        sps = config['ANNOTATOR_SPECIES'],
        homer="chisq/motif/{lib}/{sample_label}.peaks.normed.compressed.filtered.annotate.homer",
use rule peak_summary_normalized from peak_anno as chi_summary with:
    input:
        peak="chisq/{lib}/{sample_label}.peaks.normed.compressed.annotate.bed"
    output:
        summary_csv = "chisq/{lib}/{sample_label}.peaks.summary"