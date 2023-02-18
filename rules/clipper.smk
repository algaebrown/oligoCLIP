SCRIPT_PATH=config['SCRIPT_PATH']
rule clipper:
    input:
        bam="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
        bai="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam.bai"
    output:
        peak="output/CLIPper/{libname}.{sample_label}.peaks.bed"
        
    params:
        run_time="16:00:00",
        species=config['SPECIES'],
        datadir = config['DATADIR'],
        error_out_file = "error_files/clipper.{libname}.{sample_label}",
        out_file="stdout/clipper.{libname}.{sample_label}",
        cores = "10",
    benchmark: "benchmarks/clipper/clipper.{libname}.{sample_label}"
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
    benchmark: "benchmarks/clipper/count_mapped.{something}"
    shell:
        """
        module load samtools;
        samtools view -cF 4 {input} > {output};
        """
    

############# NORMALIZE TO IgG ######################
rule norm_peaks_to_another_library_inside:
    input:
        subsample_bam_ip="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
        subsample_bam_in="{libname}/bams/{bg_sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
        nread_ip="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.readnum.txt",
        nread_in="{libname}/bams/{bg_sample_label}.rmDup.Aligned.sortedByCoord.out.readnum.txt",
        peak="output/CLIPper/{libname}.{sample_label}.peaks.bed"
    output:
        norm_peak="output/CLIPper.{bg_sample_label}/{libname}.{sample_label}.peaks.normed.bed"
    params:
        run_time="4:00:00",
        error_out_file = "error_files/norm_peaks.{libname}.{sample_label}",
        out_file="stdout/norm_peaks.{libname}.{sample_label}",
        cores = "4",
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
        libname = ["{libname}"],
        sample_label = list(set(config['rbps'])-set([wildcards.clip_sample_label])-set([config['AS_INPUT']]))
        ),
    output:
        combined_bam = temp("output/CLIPper/CC_bams/{libname}.{clip_sample_label}.bam"),
        combined_sorted_bam = "output/CLIPper/CC_bams/{libname}.{clip_sample_label}.sorted.bam",
        combined_sorted_index = "output/CLIPper/CC_bams/{libname}.{clip_sample_label}.sorted.bam.bai"
    params:
        run_time="2:00:00",
        error_out_file = "error_files/combine_bam_as_CC.{libname}.{clip_sample_label}",
        out_file = "stdout/combine_bam_as_CC.{libname}.{clip_sample_label}",
        cores = "1",
    benchmark: "benchmarks/clipper/concat_bams.{libname}.{clip_sample_label}"
    shell:
        """
        module load samtools
        samtools merge {output.combined_bam} {input.other_bams}
        samtools sort {output.combined_bam} | samtools view -Sb - > {output.combined_sorted_bam}
        samtools index {output.combined_sorted_bam}
        """

use rule norm_peaks_to_another_library_inside as norm_peaks_to_cc with:
    input:
        subsample_bam_ip="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
        subsample_bam_in="output/CLIPper/CC_bams/{libname}.{sample_label}.sorted.bam",
        nread_ip="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.readnum.txt",
        nread_in="output/CLIPper/CC_bams/{libname}.{sample_label}.sorted.readnum.txt",
        peak="output/CLIPper/{libname}.{sample_label}.peaks.bed"
    output:
        norm_peak="output/CLIPper-CC/{libname}.{sample_label}.peaks.normed.bed"
    benchmark: "benchmarks/clipper/norm_peaks_cc.{libname}.{sample_label}"

rule compress_peak:
    input:
        norm_peak="{something}.normed.bed"
    output:
        compress_peak="{something}.normed.compressed.bed"
    params:
        run_time="3:00:00",
        error_out_file = "error_files/compress_peak.{something}",
        out_file="stdout/compress_peak.{something}",
        cores = "1",
    benchmark: "benchmarks/clipper/compress_peaks.{something}"
    shell:
        """
        module load eclip/0.7.0;
        perl {SCRIPT_PATH}/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.pl \
        {input.norm_peak} \
        {output.compress_peak}
        """


