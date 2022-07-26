# This is the file that annotate and run motif analysis for peaks.
SCRIPT_PATH = config['SCRIPT_PATH']

rule annotate:
    input:
        peak="output/{sample_label}.peaks.normed.compressed.bed"
    output:
        "output/{sample_label}.peaks.normed.compressed.annotate.bed"
    params:
        run_time=3,
        error_out_file = "error_files/annotate",
        cores = "1",
        gtf_db = config['GTF'],
        sps = config['ANNOTATOR_SPECIES']
    shell:
        """
        module load annotator
        annotator --input {input.peak} --output {output} --gtfdb {params.gtf_db} --species {params.sps}
        """


rule calc_partition_nuc:
    input:
        partition = "output/{sample_label}.peaks.normed.compressed.bed",
        genome = config['GENOMEFA']
    output:
        nuc = "output/{sample_label}.peaks.normed.compressed.nuc",
        gc = "output/{sample_label}.peaks.normed.compressed.gc",
    params:
        error_out_file = "error_files/partition_nuc",
        run_time = "1",
        memory = "1000",
        job_name = "calc_partition_nuc",
        cores=1
    benchmark: "benchmarks/partition_nuc.{sample_label}.txt"
    shell:
        """ 
        module load bedtools;
        bedtools nuc -s -fi {input.genome} -bed {input.partition} > {output.nuc};
        cut -d $'\t' -f 1,2,3,8 {output.nuc} > {output.gc}
        """        


rule filter_peak_fc_pval:
    input:
        "output/{sample_label}.peaks.normed.compressed.bed"
    output:
        "output/{sample_label}.peaks.normed.compressed.filtered.bed"
    conda:
        "envs/metadensity.yaml"
    params:
        run_time=1,
        error_out_file = "error_files/annotate",
        cores = "1",
    shell:
        """
        python /home/hsher/project_another/FMRP_UABP2L/renormalize_snake/filter_peak.py {input} {output}
        """

rule motif_analysis:
    input:
        peak="output/{sample_label}.peaks.normed.compressed.filtered.bed"
    output:
        pickle="output/motif/{sample_label}.peaks.normed.compressed.filtered.annotate.pickle",
        svg="output/motif/{sample_label}.peaks.normed.compressed.filtered.annotate.svg",
    params:
        run_time=3,
        error_out_file = "error_files/annotate",
        cores = "1",
        fa = config['GENOMEFA'],
        sps = config['ANNOTATOR_SPECIES'],
        homer="output/motif/{sample_label}.peaks.normed.compressed.filtered.annotate.homer",
    shell:
        """
        module load eclipanalysis
        analyze_motifs --peaks {input.peak} --k 6 \
            --out_pickle_file {output.pickle} \
            --out_homer_dir {params.homer} \
            --out_file {output.svg} \
            --species {params.sps} \
            --genome_fasta {params.fa}

        """
########### CLIPper output
 #################
rule annotate_clipper:
    input:
        peak="output/CLIPper/{sample_label}.peaks.bed"
    output:
        "output/CLIPper/{sample_label}.peaks.annotate.bed"
    params:
        run_time=3,
        error_out_file = "error_files/annotate",
        cores = "1",
        gtf_db = config['GTF'],
        sps = config['ANNOTATOR_SPECIES']
    shell:
        """
        module load annotator
        annotator --input {input.peak} --output {output} --gtfdb {params.gtf_db} --species {params.sps}
        """

rule gc_clipper:
    input:
        partition = "output/CLIPper/{sample_label}.peaks.bed",
        genome = config['GENOMEFA']
    output:
        nuc = "output/CLIPper/{sample_label}.peaks.nuc",
        gc = "output/CLIPper/{sample_label}.peaks.gc",
    params:
        error_out_file = "error_files/partition_nuc",
        run_time = "1",
        memory = "1000",
        cores=1,
        job_name = "calc_partition_nuc"
    benchmark: "benchmarks/partition_nuc.clipper.{sample_label}.txt"
    shell:
        """ 
        module load bedtools;
        bedtools nuc -s -fi {input.genome} -bed {input.partition} > {output.nuc};
        cut -d $'\t' -f 1,2,3,8 {output.nuc} > {output.gc}
        """ 

rule filter_clipper:
    input:
        "output/CLIPper/{sample_label}.peaks.bed"
    output:
        "output/CLIPper/{sample_label}.peaks.filtered.bed"
    params:
        run_time=1,
        error_out_file = "error_files/annotate",
        cores = "1",
        thres=config['CLIPper_pvalthes']
    shell:
        """
        awk '{{ if ($5 < {params.thres}) print }}' {input}  > {output}
        """

rule motif_analysis_clipper:
    input:
        peak="output/CLIPper/{sample_label}.peaks.filtered.bed"
    output:
        pickle="output/CLIPper/{sample_label}.peaks.filtered.pickle",
        svg="output/CLIPper/{sample_label}.peaks.filtered.svg",
    params:
        run_time=12,
        error_out_file = "error_files/annotate",
        cores = "1",
        fa = config['GENOMEFA'],
        sps = config['ANNOTATOR_SPECIES'],
        homer="output/CLIPper/{sample_label}.peaks.filtered.homer",
    shell:
        """
        module load eclipanalysis
        analyze_motifs --peaks {input.peak} --k 6 \
            --out_pickle_file {output.pickle} \
            --out_homer_dir {params.homer} \
            --out_file {output.svg} \
            --species {params.sps} \
            --genome_fasta {params.fa}
        """

rule peak_summary_unnormalized:
    input:
        peak="output/CLIPper/{sample_label}.peaks.annotate.bed"
    output:
        summary_csv = "summary/CLIPper/{sample_label}.peaks.summary"
    params:
        run_time=1,
        error_out_file = "error_files/annotate",
        cores = "1",
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/summarize_peak.py --bed {input.peak}  \
                                        --outfile {output.summary_csv}
        """
rule peak_summary_normalized:
    input:
        peak="output/{sample_label}.peaks.normed.compressed.annotate.bed"
    output:
        summary_csv = "summary/normalized/{sample_label}.peaks.summary"
    params:
        run_time=1,
        error_out_file = "error_files/annotate",
        cores = "1",
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/summarize_peak.py --bed {input.peak}  \
                                        --is_normalized \
                                        --outfile {output.summary_csv}
        """
