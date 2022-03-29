import pandas as pd
import os
import sys
import glob
from snakeDros_config import *
#snakemake -s snakeDros.py -j 16 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --directory=/home/hsher/scratch/katie_drosphila
 

include: "snakeDros_config.py"

# if not os.path.exists(MANIFEST): make_meta(MANIFEST)
# if EXE_DIR not in sys.path: os.environ["PATH"] = EXE_DIR + os.pathsep + os.environ["PATH"]
# if CLIP_TOOLS not in sys.path: os.environ["PATH"] = CLIP_TOOLS + os.pathsep + os.environ["PATH"]


manifest = pd.read_csv(MANIFEST, index_col = False)

sample_labels = manifest.Sample.tolist()
print(manifest.iloc[0])
manifest.loc[manifest.Sample == sample_labels[0]]["RNA_adaptor"].values[0]

print(os.path.join(ADAPTOR_PATH, 
        manifest.loc[manifest.Sample == sample_labels[0], "RNA_adaptor"].values[0]+'_adapters.fasta'))

rule all:
    input:
        expand("bams/repeat/{sample_label}.Aligned.out.sort.bam.bai", sample_label = sample_labels)+
        expand("fastqc/{sample_label}.umi.fqTrTr_fastqc.html", sample_label = sample_labels)+
        expand("bams/genome/{sample_label}.genome-mappedSoSo.rmDupSo.Aligned.out.bam.bai", sample_label = sample_labels)
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





rule extract_umi:
    input:
        fq_raw = lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["fastq_1"].values[0]),
    output:
        fq_umi = "fastqs/umi/{sample_label}.umi.fq.gz",
        metrics = "{sample_label}.umi.metrics"
    params:
        error_out_file = "error_files/{sample_label}_extract_sra",
        run_time = "1:45:00",
        cores = "4",
        memory = "10000",
        job_name = "extract_umi",
        umi_pattern = lambda wildcards: 'N'*(manifest.loc[manifest.Sample == wildcards.sample_label]["random_mer"].values[0]),
    benchmark: "benchmarks/umi/{sample_label}.extract.txt"
    shell:
        """
        module load eclip;
        umi_tools extract \
            --random-seed 1 \
            --bc-pattern {params.umi_pattern} \
            --stdin {input.fq_raw} \
            --stdout {output.fq_umi} \
            --log {output.metrics}
        """

#TODO: make "find_InvRNAadaptor"
#TODO: try if 1 step is good enough
#TODO: try if we can feed the entire barcode without tiling

def get_full_adapter_path(sample_label):
    return os.path.join(ADAPTOR_PATH, manifest.loc[manifest.Sample == sample_label, "RNA_adaptor"].values[0]+'_adapters.fasta')

rule cutadapt_round_one:
    input:
        fq_umi = "fastqs/umi/{sample_label}.umi.fq.gz",
    output:
        fq_trimmed="fastqs/umi/{sample_label}.umi.fqTr.gz",
        metrics = "fastqs/umi/{sample_label}.umi.r1.fqTr.metrics"
    params:
        InvRNA=lambda wildcards: get_full_adapter_path(wildcards.sample_label),
        run_time = "08:04:00",
        cores="4"
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
        fq_trimmed="fastqs/umi/{sample_label}.umi.fqTr.gz",
    output:
        fq_trimmed_twice="fastqs/umi/{sample_label}.umi.fqTrTr.gz",
        metrics = "fastqs/umi/{sample_label}.umi.r1.fqTrTr.metrics"
    params:
        InvRNA=lambda wildcards: get_full_adapter_path(wildcards.sample_label),
        run_time = "08:04:00",
        cores="4"
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
        fq_trimmed_twice="fastqs/umi/{sample_label}.umi.fqTrTr.gz",
    output:
        fq="fastqs/umi/{sample_label}.umi.fqTrTr",
        fq_gz="fastqs/umi/{sample_label}.umi.fqTrTr.sorted.fq.gz",
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
        "fastqs/umi/{sample_label}.umi.fqTrTr.gz"
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
rule align_reads_to_Drosophila:
    input:
        fq_1 = "fastqs/umi/{sample_label}.umi.fqTrTr.sorted.fq.gz",
    output:
        ubam = "bams/dros/{sample_label}_dros_mapped.Aligned.out.bam",
        unmqpped= "bams/dros/{sample_label}_dros_mapped.Unmapped.out.mate1",
        log= "bams/dros/{sample_label}_dros_mapped.Log.final.out",
    params:
        error_out_file = "error_files/{sample_label}_align_dros_reads",
        run_time = "08:40:00",
        cores = "4",
        memory = "10000",
        job_name = "align_reads",
        star_sjdb = STAR_DROS,
        outprefix = "bams/dros/{sample_label}_dros_mapped.",
    benchmark: "benchmarks/align/{sample_label}.align_dros_reads.txt"
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

rule align_reads_to_REPEAT:
    input:
        fq_1 = "bams/dros/{sample_label}_dros_mapped.Unmapped.out.mate1"
    output:
        ubam = "bams/repeat/{sample_label}.Aligned.out.bam",
        unmqpped= "bams/repeat/{sample_label}.Unmapped.out.mate1",
        log= "bams/repeat/{sample_label}.Log.final.out",
    params:
        error_out_file = "error_files/{sample_label}_align_reads",
        run_time = "08:40:00",
        cores = "4",
        memory = "10000",
        job_name = "align_reads",
        star_sjdb = STAR_REP,
        outprefix = "bams/repeat/{sample_label}.",
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
            --runMode alignReads \
            --runThreadN 8
        """
rule align_to_GENOME:
    input:
        fq= "bams/repeat/{sample_label}.Unmapped.out.mate1",
    output:
        ubam = "bams/genome/{sample_label}.genome-mapped.Aligned.out.bam",
        unmqpped= "bams/genome/{sample_label}.genome-mapped.Unmapped.out.mate1",
        log= "bams/genome/{sample_label}.genome-mapped.Log.final.out",
    params:
        error_out_file = "error_files/{sample_label}_align_reads_genome",
        run_time = "08:40:00",
        cores = "4",
        memory = "10000",
        job_name = "align_reads",
        star_sjdb = STAR_DIR,
        outprefix = "bams/genome/{sample_label}.genome-mapped.",
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
        bam="bams/genome/{sample_label}.genome-mapped.Aligned.out.bam",
    output:
        sort_once = "bams/genome/{sample_label}.genome-mappedSo.Aligned.out.bam",
        sort_twice= "bams/genome/{sample_label}.genome-mappedSoSo.Aligned.out.bam",
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
        """

rule index_bams:
    input:
        ubam = "bams/repeat/{sample_label}.Aligned.out.bam"
    output:
        sbam = "bams/repeat/{sample_label}.Aligned.out.sort.bam",
        ibam = "bams/repeat/{sample_label}.Aligned.out.sort.bam.bai"
    params:
        error_out_file = "error_files/{sample_label}_index_bams",
        run_time = "01:40:00",
        cores = "4",
        memory = "1000",
        job_name = "index_bam"
    benchmark: "benchmarks/align/{sample_label}.index_bam.txt"
    shell:
        "module load samtools;"
        "samtools sort -T {wildcards.sample_label} -@ 4 {input.ubam} > {output.sbam};"
        "samtools index {output.sbam};"

rule index_bams2:
    input:
        bam = "bams/genome/{sample_label}.genome-mappedSoSo.Aligned.out.bam"
    output:
        ibam = "bams/genome/{sample_label}.genome-mappedSoSo.Aligned.out.bam.bai"
    params:
        error_out_file = "error_files/{sample_label}_index_bams",
        run_time = "01:40:00",
        cores = "4",
        memory = "1000",
        job_name = "index_bam"
    benchmark: "benchmarks/align/{sample_label}.index_bam.txt"
    shell:
        "module load samtools;"
        "samtools index {input.bam};"


rule umi_dedup:
    input:
        bam="bams/genome/{sample_label}.genome-mappedSoSo.Aligned.out.bam",
        bai="bams/genome/{sample_label}.genome-mappedSoSo.Aligned.out.bam.bai",
    output:
        bam_dedup="bams/genome/{sample_label}.genome-mappedSoSo.rmDup.Aligned.out.bam"
    params:
        error_out_file = "error_files/{sample_label}_sort_bam",
        run_time = "03:40:00",
        cores = "4",
        memory = "10000",
        job_name = "sortbam",
        prefix='bams/genome/{sample_label}.genome-mappedSoSo'
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
        bam = "bams/genome/{sample_label}.genome-mappedSoSo.rmDup.Aligned.out.bam"
    output:
        sbam="bams/genome/{sample_label}.genome-mappedSoSo.rmDupSo.Aligned.out.bam",
        bai = "bams/genome/{sample_label}.genome-mappedSoSo.rmDupSo.Aligned.out.bam.bai"
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
