import pandas as pd
import os
import sys
import glob
from snakeCLIP_config import *
#snakemake -s snakeCLIP.py -j 10 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q glean" --directory /home/hsher/scratch/ABC

include: "snakeCLIP_config.py"
rbps = pd.read_csv(barcode, header = None, sep = '\t', names = ['RBP', 'barcode'])['RBP'].tolist()


rule all:
    input:
        expand("fastqc/{sample_label}.umi.fqTrTr_fastqc.html", sample_label = rbps)+
        expand("bams/{sample_label}.genome-mappedSoSo.rmDupSo.Aligned.out.bam.bai", sample_label = rbps),
        "whitelist.txt",
        "merged.csv"
    output:
        "snakeLeaf.txt"
    params:
        error_out_file = "error_files/all",
        run_time = "00:04:00",
        cores = "1",
        memory = "20",
        job_name = "all"
    shell:
        "echo $(date)  > {output};"
        "echo created by Evan Boyle and the Yeo lab >> {output}"

def fastqc_name(full_path):
    return os.path.basename(full_path).replace('.fastqc.gz', '_fastqc.html')

# rule fastqc_initial:
#     input:
#         lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["fastq_1"].values[0])
#     output:
#         "fastqc/{sample_label}_fastqc.html"
#     threads: 2
#     params:
#         outdir="fastqc/{sample_label}",
#         original_name = lambda wildcards: fastqc_name(manifest.loc[manifest.Sample == wildcards.sample_label]["fastq_1"].values[0]),
#         cores="1",
#         run_time = "00:04:00",
#     resources:
#         runtime="1:00:00",
        
#     shell:
#         """
#         module load fastqc;
#         fastqc {input} --extract --outdir {params.outdir} -t {threads}
#         mv fastqc/{params} {output}
#         """


rule demultiplex:
    input:
        fq_raw = FASTQ,
        barcode_tsv= barcode
    output:
        expand("{sample_label}.fastq", sample_label = rbps)
    params:
        outdir="",
        cores="1",
        run_time = "03:00:00",
    shell:
        """
        module load fastx_toolkit
        zcat {input.fq_raw} | fastx_barcode_splitter.pl --bcfile {input.barcode_tsv} --prefix "" --suffix ".fastq" --bol
        """

rule find_barcode:
    input:
        fq_raw = FASTQ
    output:
        "whitelist.txt"
    params:
        error_out_file = "error_files/findbarcode",
        run_time = "1:45:00",
        cores = "4",
        memory = "10000",
        job_name = "extract_umi",
        umi_pattern = umi_pattern,
        n_cell = cell_number
    shell:
        """
        module load eclip;
        umi_tools whitelist --stdin {input.fq_raw} --bc-pattern {params.umi_pattern} --set-cell-number={params.n_cell} > {output}
        """

rule extract_umi:
    input:
        barcode = "whitelist.txt",
        fq_raw = "{sample_label}.fastq",
    output:
        fq_umi = "fastqs/{sample_label}.umi.fq.gz",
        metrics = "fastqs/{sample_label}.umi.metrics"
    params:
        error_out_file = "error_files/{sample_label}_extract_sra",
        run_time = "1:45:00",
        cores = "4",
        memory = "10000",
        job_name = "extract_umi",
        umi_pattern = umi_pattern
    benchmark: "benchmarks/umi/{sample_label}.extract.txt"
    shell:
        """
        module load eclip;
        umi_tools extract \
            --random-seed 1 \
            --bc-pattern {params.umi_pattern} \
            --stdin {input.fq_raw} \
            --stdout {output.fq_umi} \
            --log {output.metrics} \
            --whitelist={input.barcode}
        """

# TODO: generate warning when some barcode is non-exsitence
rule combine_barcode:
    input:
        whitelist = "whitelist.txt",
        provided_barcode = barcode
    output:
        "merged.csv"
    params:
        error_out_file = "error_files/merge",
        run_time = "0:30:00",
        cores = "1",
        python = PYTHON3_PATH
    shell:
        """
        {params.python} scripts/barcode_join.py {input.provided_barcode} {input.whitelist} {output}
        """


#TODO: make "find_InvRNAadaptor"
#TODO: try if 1 step is good enough
#TODO: try if we can feed the entire barcode without tiling

def get_full_adapter_path(adaptor):
    return os.path.join(ADAPTOR_PATH, adaptor+'_adapters.fasta')

rule cutadapt_round_one:
    input:
        fq_umi = "fastqs/{sample_label}.umi.fq.gz",
    output:
        fq_trimmed="fastqs/{sample_label}.umi.fqTr.gz",
        metrics = "fastqs/{sample_label}.umi.r1.fqTr.metrics"
    params:
        InvRNA=lambda wildcards: get_full_adapter_path(adaptor),
        run_time = "02:04:00",
        cores="1"
    benchmark: "benchmarks/cutadapt/{sample_label}.extract.txt"
    shell:
        """
        module load eclip;
        cutadapt -O 1 \
            -f fastq \
            --match-read-wildcards \
            --times 1 \
            -e 0.1 \
            --quality-cutoff 6 \
            -m 18 \
            -o {output.fq_trimmed} \
            -a file:{params.InvRNA} \
            {input.fq_umi} > {output.metrics}
        """

rule cutadapt_round_two:
    input:
        fq_trimmed="fastqs/{sample_label}.umi.fqTr.gz",
    output:
        fq_trimmed_twice="fastqs/{sample_label}.umi.fqTrTr.gz",
        metrics = "fastqs/{sample_label}.umi.r1.fqTrTr.metrics"
    params:
        InvRNA=lambda wildcards: get_full_adapter_path(adaptor),
        run_time = "02:04:00",
        cores="1"
    benchmark: "benchmarks/cutadapt/{sample_label}.extract_round2.txt"
    shell:
        """
        module load eclip;
        cutadapt -O 5 \
            -f fastq \
            --match-read-wildcards \
            --times 1 \
            -e 0.1 \
            --quality-cutoff 6 \
            -m 18 \
            -o {output.fq_trimmed_twice}\
            -a file:{params.InvRNA} \
            {input.fq_trimmed} > {output.metrics}
        """
rule sort_and_gzip_fastq:
    input:
        fq_trimmed_twice="fastqs/{sample_label}.umi.fqTrTr.gz",
    output:
        fq="fastqs/{sample_label}.umi.fqTrTr",
        fq_gz="fastqs/{sample_label}.umi.fqTrTr.sorted.fq.gz",
    params:
        run_time = "05:04:00",
        cores="1"
    shell:
        """
        module load eclip;
        zcat {input.fq_trimmed_twice} > {output.fq}
        fastq-sort --id {output.fq} | gzip > {output.fq_gz}
        """
# TODO, CHECK IF THE TRIMMING IS SUCCESSFUL, CHECK CROSS CONTAMINATION
rule fastqc_post_trim:
    input:
        "fastqs/{sample_label}.umi.fqTrTr.gz"
    output:
        "fastqc/{sample_label}.umi.fqTrTr_fastqc.html"
    threads: 2
    params:
        outdir="fastqc",
        run_time = "01:09:00",
        cores="1"
    resources:
        runtime="1:00:00",
        cores="1"
    shell:
        """
        module load fastqc;
        fastqc {input} --extract --outdir {params.outdir} -t {threads}
        """

rule align_reads_to_REPEAT:
    input:
        fq_1 = "fastqs/{sample_label}.umi.fqTrTr.sorted.fq.gz",
    output:
        ubam = "bams/{sample_label}.Aligned.out.bam",
        unmqpped= "bams/{sample_label}.Unmapped.out.mate1",
        log= "bams/{sample_label}.Log.final.out",
    params:
        error_out_file = "error_files/{sample_label}_align_reads",
        run_time = "03:40:00",
        cores = "4",
        memory = "10000",
        job_name = "align_reads",
        star_sjdb = STAR_REP,
        outprefix = "bams/{sample_label}.",
    benchmark: "benchmarks/align/{sample_label}.align_reads.txt"
    shell:
        """
        module load star ;
        STAR \
            --alignEndsType EndToEnd \
            --genomeDir {params.star_sjdb} \
            --genomeLoad NoSharedMemory \
            --outBAMcompression 10 \
            --outFileNamePrefix {params.outprefix} \
            --outFilterMultimapNmax 1 \
            --outFilterMultimapScoreRange 1 \
            --outFilterScoreMin 10 \
            --outFilterType BySJout \
            --outReadsUnmapped Fastx \
            --outSAMattrRGline ID:foo \
            --outSAMattributes All \
            --outSAMmode Full \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --outStd Log \
            --readFilesIn {input.fq_1} \
            --readFilesCommand zcat \
            --runMode alignReads \
            --runThreadN 8
        """
rule align_to_GENOME:
    input:
        fq= "bams/{sample_label}.Unmapped.out.mate1",
    output:
        ubam = "bams/{sample_label}.genome-mapped.Aligned.out.bam",
        unmqpped= "bams/{sample_label}.genome-mapped.Unmapped.out.mate1",
        log= "bams/{sample_label}.genome-mapped.Log.final.out",
    params:
        error_out_file = "error_files/{sample_label}_align_reads_genome",
        run_time = "03:40:00",
        cores = "4",
        memory = "10000",
        job_name = "align_reads",
        star_sjdb = STAR_DIR,
        outprefix = "bams/{sample_label}.genome-mapped.",
    benchmark: "benchmarks/align/{sample_label}.align_reads.txt"
    shell:
        """
        module load star ;
        STAR \
        --alignEndsType EndToEnd \
        --genomeDir {params.star_sjdb} \
        --genomeLoad NoSharedMemory \
        --outBAMcompression 10 \
        --outFileNamePrefix {params.outprefix} \
        --outFilterMultimapNmax 1 \
        --outFilterMultimapScoreRange 1 \
        --outFilterScoreMin 10 \
        --outFilterType BySJout \
        --outReadsUnmapped Fastx \
        --outSAMattrRGline ID:foo \
        --outSAMattributes All \
        --outSAMmode Full \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outStd Log \
        --readFilesIn {input.fq} \
        --runMode alignReads \
        --runThreadN 8

        """
rule sort_bams:
    input:
        bam="bams/{sample_label}.genome-mapped.Aligned.out.bam",
    output:
        sort_once = "bams/{sample_label}.genome-mappedSo.Aligned.out.bam",
        sort_twice= "bams/{sample_label}.genome-mappedSoSo.Aligned.out.bam",
        bai = "bams/{sample_label}.genome-mappedSoSo.Aligned.out.bam.bai"
    params:
        error_out_file = "error_files/{sample_label}_sort_bam",
        run_time = "03:40:00",
        cores = "4",
        memory = "10000",
        job_name = "sortbam",
    shell:
        """
        module load samtools;
        samtools sort -o {output.sort_once} {input.bam};
        samtools sort -o {output.sort_twice} {output.sort_once};
        samtools index {output.sort_twice}
        """


rule umi_dedup:
    input:
        bam="bams/{sample_label}.genome-mappedSoSo.Aligned.out.bam",
        bai="bams/{sample_label}.genome-mappedSoSo.Aligned.out.bam.bai"
    output:
        bam_dedup="bams/{sample_label}.genome-mappedSoSo.rmDup.Aligned.out.bam"
    params:
        error_out_file = "error_files/{sample_label}_sort_bam",
        run_time = "03:40:00",
        cores = "4",
        memory = "10000",
        job_name = "sortbam",
        prefix='bams/{sample_label}.genome-mappedSoSo'
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
rule index_genome_bams:
    input:
        bam = "bams/{sample_label}.genome-mappedSoSo.rmDup.Aligned.out.bam"
    output:
        sbam="bams/{sample_label}.genome-mappedSoSo.rmDupSo.Aligned.out.bam",
        bai = "bams/{sample_label}.genome-mappedSoSo.rmDupSo.Aligned.out.bam.bai"
    params:
        error_out_file = "error_files/{sample_label}_index_bams",
        run_time = "01:40:00",
        cores = "4",
        memory = "1000",
        job_name = "index_bam"
    benchmark: "benchmarks/align/{sample_label}.index_bam.txt"
    shell:
        "module load samtools;"
        "samtools sort -o {output.sbam} {input.bam} ;"
        "samtools index {output.sbam};"
