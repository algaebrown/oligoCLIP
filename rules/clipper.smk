locals().update(config)
# rule make_clipper_annotation:
#     input:
#         config['GTF'] # gtf.db
#     output:

#     params:

#     shell:
#         """
#         module load annotator/0.99.0
#         create_region_bedfiles \
#             --db_file {input} \  # database file downloaded from above
#             --species mm10 \  # sets the gff/gtf nomenclature (essentially whether it's gencode or wormbase gtf format)
#             --cds_out outputs/mm10_vM10_cds.bed \  # output cds region
#             --proxintron_out outputs/mm10_vM10_proxintrons.bed \  # output proximal intron regions (500nt from exons)
#             --distintron_out outputs/mm10_vM10_distintrons.bed \  # output distal intron regions (> 500nt from exons)
#             --utr5_out outputs/mm10_vM10_five_prime_utrs.bed \  # output 5'UTR regions
#             --utr3_out outputs/mm10_vM10_three_prime_utrs.bed  # output 3' UTR regions
#         """

rule clipper:
    input:
        bam="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
        bai="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam.bai"
    output:
        peak="CLIPper/{libname}.{sample_label}.peaks.bed"
    params:
        run_time="16:00:00",
        species=config['SPECIES'],
        datadir = config['DATADIR'],
        error_out_file = "error_files/clipper.{libname}.{sample_label}",
        out_file="stdout/clipper.{libname}.{sample_label}",
        cores = "10",
        memory = 40000,
    benchmark: "benchmarks/clipper/clipper.{libname}.{sample_label}"
    container:
        "docker://brianyee/clipper:charlene_move_data"
    shell:
        """
        module load clipper/charlene_move_data_branch
        clipper -b {input.bam} -s {params.species} -o {output.peak} --processors=16 --datadir {params.datadir}
        """

rule count_mapped_reads:
    input:
        "{something}.bam",
    output:
        "{something}.readnum.txt",
    params:
        run_time="00:30:00",
        error_out_file = "error_files/countread.{something}",
        out_file = "stdout/countread.{something}",
        cores = "1",
        memory = 1000,
    benchmark: "benchmarks/clipper/count_mapped.{something}"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools view -cF 4 {input} > {output};
        """
    

############# NORMALIZE TO Internal IgG ######################
rule norm_peaks_to_another_library_inside:
    input:
        subsample_bam_ip="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
        subsample_bam_in="{libname}/bams/{bg_sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
        nread_ip="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.readnum.txt",
        nread_in="{libname}/bams/{bg_sample_label}.rmDup.Aligned.sortedByCoord.out.readnum.txt",
        peak="CLIPper/{libname}.{sample_label}.peaks.bed"
    output:
        norm_peak=temp("CLIPper.{bg_sample_label}/{libname}.{sample_label}.peaks.normed.bed")
    params:
        run_time="4:00:00",
        error_out_file = "error_files/norm_peaks.{libname}.{sample_label}",
        out_file="stdout/norm_peaks.{libname}.{sample_label}",
        cores = "4",
        memory = 10000,
    benchmark: "benchmarks/clipper/norm_peaks.{libname}.{sample_label}.{bg_sample_label}"
    shell:
        """
        module load eclip/0.7.0;
        perl {SCRIPT_PATH}/overlap_peakfi_with_bam.pl \
            {input.subsample_bam_ip} \
            {input.subsample_bam_in} \
            {input.peak} \
            {input.nread_ip} \
            {input.nread_in} \
            {output.norm_peak}
        """

############# COMPLEMENTARY CONTROL ######################
rule concat_other_bams_as_complementary_control:
    input:
        other_bams = lambda wildcards: expand("{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
            libname = [wildcards.libname],
            sample_label = list(set(config['rbps'])-set([wildcards.clip_sample_label])-set(config['AS_INPUT']))
        ),
    output:
        combined_bam = temp("CLIPper/CC_bams/{libname}.{clip_sample_label}.bam"),
        combined_sorted_bam = temp("CLIPper/CC_bams/{libname}.{clip_sample_label}.sorted.bam"),
        combined_sorted_index = temp("CLIPper/CC_bams/{libname}.{clip_sample_label}.sorted.bam.bai")
    params:
        run_time="2:00:00",
        error_out_file = "error_files/combine_bam_as_CC.{libname}.{clip_sample_label}",
        out_file = "stdout/combine_bam_as_CC.{libname}.{clip_sample_label}",
        cores = "1",
        memory = 40000,
    benchmark: "benchmarks/clipper/concat_bams.{libname}.{clip_sample_label}"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools merge {output.combined_bam} {input.other_bams}
        samtools sort {output.combined_bam} | samtools view -Sb - > {output.combined_sorted_bam}
        samtools index {output.combined_sorted_bam}
        """

use rule norm_peaks_to_another_library_inside as norm_peaks_to_cc with:
    input:
        subsample_bam_ip="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
        subsample_bam_in="CLIPper/CC_bams/{libname}.{sample_label}.sorted.bam",
        nread_ip="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.readnum.txt",
        nread_in="CLIPper/CC_bams/{libname}.{sample_label}.sorted.readnum.txt",
        peak="CLIPper/{libname}.{sample_label}.peaks.bed",
    output:
        norm_peak=temp("CLIPper_CC/{libname}.{sample_label}.peaks.normed.bed")
    benchmark: "benchmarks/clipper/norm_peaks_cc.{libname}.{sample_label}"

############ external normalization ##################
use rule count_mapped_reads as count_mapped_reads_external with:
    input:
        lambda wildcards: config['external_bam'][wildcards.external_label]['file']
    output:
        "CLIPper/external_readnum/{external_label}.readnum.txt"
    params:
        run_time="00:30:00",
        error_out_file = "error_files/countread.{external_label}",
        out_file = "stdout/countread.{external_label}",
        cores = "1",
        memory = 1000,
    benchmark: "benchmarks/clipper/count_mapped.{external_label}"
    


use rule norm_peaks_to_another_library_inside as norm_peaks_to_external with:
    input:
        subsample_bam_ip="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
        subsample_bam_in=lambda wildcards: ancient(config['external_bam'][wildcards.external_label]['file']),
        nread_ip="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.readnum.txt",
        nread_in="CLIPper/external_readnum/{external_label}.readnum.txt",
        peak="CLIPper/{libname}.{sample_label}.peaks.bed"
    output:
        norm_peak="CLIPper-{external_label}/{libname}.{sample_label}.peaks.normed.bed"
    benchmark: "benchmarks/clipper/norm_peaks_{external_label}.{libname}.{sample_label}"

############ commonly used ##################

rule compress_peak:
    input:
        norm_peak="{something}.normed.bed"
    output:
        compress_peak=temp("{something}.normed.compressed.bed")
    params:
        run_time="3:00:00",
        error_out_file = "error_files/compress_peak.{something}",
        out_file="stdout/compress_peak.{something}",
        cores = "1",
        memory = 40000,
    benchmark: "benchmarks/clipper/compress_peaks.{something}"
    shell:
        """
        module load eclip/0.7.0;
        perl {SCRIPT_PATH}/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.pl \
        {input.norm_peak} \
        {output.compress_peak}
        """


