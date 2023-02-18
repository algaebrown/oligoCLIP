rule annotate:
    input:
        peak="{something}.normed.compressed.bed"
    output:
        "{something}.normed.compressed.annotate.bed"
    params:
        run_time="3:00:00",
        error_out_file = "error_files/annotate.{something}",
        out_file =  "stdout/annotate.{something}",
        cores = "1",
        gtf_db = config['GTF'],
        sps = config['ANNOTATOR_SPECIES']
    shell:
        """
        module load annotator
        annotator --input {input.peak} --output {output} --gtfdb {params.gtf_db} --species {params.sps}
        """
# rule calc_partition_nuc:
#     input:
#         partition = "output/{sample_label}.peaks.normed.compressed.bed",
#         genome = config['GENOMEFA']
#     output:
#         nuc = "output/{sample_label}.peaks.normed.compressed.nuc",
#         gc = "output/{sample_label}.peaks.normed.compressed.gc",
#     params:
#         error_out_file = "error_files/partition_nuc",
#         run_time = "1",
#         memory = "1000",
#         job_name = "calc_partition_nuc",
#         cores=1
#     benchmark: "benchmarks/partition_nuc.{sample_label}.txt"
#     shell:
#         """ 
#         module load bedtools;
#         bedtools nuc -s -fi {input.genome} -bed {input.partition} > {output.nuc};
#         cut -d $'\t' -f 1,2,3,8 {output.nuc} > {output.gc}
#         """

rule motif_analysis:
    input:
        peak="{something}.normed.compressed.bed"
    output:
        filtered_peak = "{something}.normed.compressed.filtered.bed",
        pickle="{something}.normed.compressed.motif.pickle",
        svg="{something}.normed.compressed.motif.svg",
    params:
        run_time="3:00:00",
        error_out_file = "error_files/motif.{something}",
        out_file =  "stdout/motif.{something}",
        cores = "1",
        fa = config['GENOMEFA'],
        sps = config['ANNOTATOR_SPECIES'],
        homer="{something}.normed.compressed.homer",
        log10p_thres=3,
        l2fc_thres = 3,
    shell:
        """
        module load eclipanalysis

        awk '{{ if ($4 > {params.log10p_thres})&($5 > {params.l2fc_thres}) print }}' {input.peak}  > {output.filtered_peak}

        analyze_motifs --peaks {input.peak} --k 6 \
            --out_pickle_file {output.pickle} \
            --out_homer_dir {params.homer} \
            --out_file {output.svg} \
            --species {params.sps} \
            --genome_fasta {params.fa}
        """