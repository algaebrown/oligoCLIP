
# this file preprocess ABC. fastq -> bam
import pandas as pd
import os
import sys
import glob
#snakemake -s snakeCLIP.py -j 12 --keep-going --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --directory /home/hsher/scratch/ABC_reprocess/ {K562_SLBP_rep1,K562_SLBP_rep2}/bams/genome/SLBP.genome-mappedSoSo.rmDupSo.Aligned.out.bam.bai
#snakemake -s snakeCLIP.py -j 12 --keep-going --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --directory /home/hsher/scratch/ABC_reprocess/ {K562_RBFOX2_rep1,K562_RBFOX2_rep2}/bams/genome/RBFOX2.genome-mappedSoSo.rmDupSo.Aligned.out.bam.bai --configfile eclipse_rbfox.yaml
# snakemake -s snakeCLIP.py -j 12 --keep-going --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --directory /home/hsher/scratch/ABC_reprocess/ --configfile config/preprocess_config/eclipse_multi.yaml -n
# snakemake -s snakeCLIP.py -j 12 --keep-going --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --configfile config/preprocess_config/eclipse_multi2.yaml -n
#snakemake -s snakeCLIP.py -j 12 --keep-going --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --directory /home/hsher/scratch/ABC_katie/ --configfile eclipse_rbfox_katie.yaml --use-conda
#snakemake -s snakeCLIP.py -j 12 --keep-going --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --directory /home/hsher/scratch/ABC_katie/ --configfile eclipse_igg_katie.yaml --use-conda
#snakemake -s snakeCLIP.py -j 12 --keep-going --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --directory /home/hsher/scratch/ABC_reprocess/ --configfile eclipse_rbfox2_singleplex_rep.yaml --use-conda -n


try:
    print('overriding with config')
    fastq_menifest = pd.read_csv(config['fastq_menifest'])
    barcode = config['barcode']

    ADAPTOR_PATH = config['ADAPTOR_PATH']
    adaptor = config['adaptor']

    # Resources
    CHROM_SIZES = config['CHROM_SIZES']
    STAR_DIR = config['STAR_DIR']
    STAR_REP= config['STAR_REP']
    
    TOOL_PATH= config['TOOL_PATH']
    umi_pattern = config['umi_pattern']

    workdir: config['WORKDIR']

    STAR_DROS = config.get('STAR_DROS')
    
except Exception as e:
    print(e)



libs = fastq_menifest['libname'].tolist()
print(fastq_menifest)
print('LIBRARY:',libs)
barcode_df = pd.read_csv(barcode)
rbps = barcode_df['rbp'].tolist()
print('RBP:',rbps)
print('STAR:DROS',STAR_DROS)

module snakeDros:
    snakefile:
        "snakeDros.py"

module QC:
    snakefile:
        "rules/QC.py"
    config:
        config

def get_output():
    output = expand("{libname}/fastqc/{sample_label}.umi.fqTrTr.rev.sorted_fastqc.html", libname = libs, sample_label = rbps
    )+expand("{libname}/bams/genome/{sample_label}.genome-mappedSoSo.rmDupSo.Aligned.out.bam.bai",  libname = libs, sample_label = rbps
    )+expand("QC/repeat_mapping_stats.csv",  libname = libs
    )+expand("QC/genome_mapping_stats.csv",  libname = libs
    )+expand('QC/{libname}/fastQC_basic_summary.csv',  libname = libs
    )+expand('QC/{libname}/fastQC_passfail.csv',  libname = libs
    )+['QC/cutadapt_log1.csv','QC/cutadapt_log2.csv']

    if STAR_DROS:
        output.append(expand('QC/dros_mapping_stats.csv',  libname = libs))
    
    return output

rule all:
    input:
        get_output()
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

rule make_good_barcode_tsv:
    input:
        barcode
    output:
        "{libname}/barcode.tsv"
    params:
        error_out_file = "error_files/all",
        run_time = "00:04:00",
        cores = "1",
        memory = "20",
        job_name = "all",
        script_path = os.path.join(TOOL_PATH, 'scripts/to_tsv.py')
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {params.script_path} {input} {output}
        """

rule extract_umi: # TODO: adaptor TRIM first
    input:
        fq_raw = lambda wildcards: fastq_menifest.loc[fastq_menifest['libname']==wildcards.libname, 'fastq'].iloc[0]
    output:
        fq_umi = "{libname}/fastqs/all.umi.fq.gz",
        metrics = "QC/{libname}.all.umi.metrics"
    params:
        error_out_file = "error_files/extract_umi",
        run_time = "3:45:00",
        cores = "4",
        memory = "10000",
        job_name = "extract_umi",
        umi_pattern = umi_pattern
    benchmark: "benchmarks/umi/extract.{libname}.txt"
    shell:
        """
        module load eclip;
        umi_tools extract \
            --random-seed 1 \
            --bc-pattern {params.umi_pattern} \
            --stdin {input.fq_raw} \
            --stdout {output.fq_umi} \
            --log {output.metrics} \
        """

def get_full_adapter_path(adaptor):
    return os.path.join(ADAPTOR_PATH, adaptor+'_adapters.fasta')

rule cutadapt_round_one:
    input:
        fq_umi = "{libname}/fastqs/all.umi.fq.gz",
    output:
        fq_trimmed="{libname}/fastqs/all.umi.fqTr.gz",
        metrics = "QC/{libname}.umi.r1.fqTr.metrics"
    params:
        InvRNA=lambda wildcards: get_full_adapter_path(adaptor),
        run_time = "12:04:00",
        cores="4"
    benchmark: "benchmarks/cutadapt/extract.{libname}.txt"
    shell:
        """
        module load cutadapt/2.8;
        cutadapt -O 1 \
            --match-read-wildcards \
            --times 1 \
            -e 0.1 \
            --quality-cutoff 6 \
            -m 23 \
            -o {output.fq_trimmed} \
            -a file:{params.InvRNA} \
            --cores=0 \
            {input.fq_umi} > {output.metrics}
        """

rule cutadapt_round_two:
    input:
        fq_trimmed="{libname}/fastqs/all.umi.fqTr.gz",
    output:
        fq_trimmed_twice="{libname}/fastqs/all.umi.fqTrTr.gz",
        metrics = "QC/{libname}.umi.r1.fqTrTr.metrics"
    params:
        InvRNA=lambda wildcards: get_full_adapter_path(adaptor),
        run_time = "16:04:00",
        cores="8"
    benchmark: "benchmarks/cutadapt/extract_round2.{libname}.txt"
    shell:
        """
        module load cutadapt/2.8;
        cutadapt -O 5 \
            --match-read-wildcards \
            --times 1 \
            -e 0.1 \
            --quality-cutoff 6 \
            -m 23 \
            -o {output.fq_trimmed_twice}\
            -a file:{params.InvRNA} \
            --cores=0 \
            {input.fq_trimmed} > {output.metrics}
        """

use rule gather_trimming_stat from QC as qc_trim1 with:
    input:
        tr1=expand("QC/{libname}.umi.r1.fqTr.metrics", libname = libs),
    output:
        tr1='QC/cutadapt_log1.csv'

use rule gather_trimming_stat from QC as qc_trim2 with:
    input:
        tr1=expand("QC/{libname}.umi.r1.fqTrTr.metrics", libname = libs)
    output:
        tr1='QC/cutadapt_log2.csv',

rule demultiplex:
    input:
        fq_raw = "{libname}/fastqs/all.umi.fqTrTr.gz",
        barcode_tsv= "{libname}/barcode.tsv",
    output:
        fq=expand("{libname}/fastqs/{sample_label}.umi.fqTrTr.fastq", libname = ["{libname}"], sample_label = rbps),
        logs = "{libname}/barcode.log"
    params:
        outdir="",
        cores="1",
        run_time = "03:00:00",
        prefix = "{libname}/fastqs/"
    shell:
        """
        module load fastx_toolkit
        zcat {input.fq_raw} | fastx_barcode_splitter.pl --bcfile {input.barcode_tsv} --prefix {params.prefix} --suffix ".umi.fqTrTr.fastq" --bol > {output.logs}
        """
rule remove_barcode_and_reverse_complement:
    input:
        "{libname}/fastqs/{sample_label}.umi.fqTrTr.fastq"
    output:
        fqrev="{libname}/fastqs/{sample_label}.umi.fqTrTr.rev.fq",
    params:
        run_time = "05:04:00",
        cores="1"
    shell:
        """
        module load fastx_toolkit
        module load eclip

        cutadapt -g NNNNN {input} | fastx_reverse_complement > {output.fqrev}
        """
rule sort_and_gzip_fastq:
    input:
        fq_trimmed_twice="{libname}/fastqs/{sample_label}.umi.fqTrTr.rev.fq",
    output:
        fq_gz="{libname}/fastqs/{sample_label}.umi.fqTrTr.rev.sorted.fq.gz"
    params:
        run_time = "05:04:00",
        cores="1"
    shell:
        """
        module load eclip;
        fastq-sort --id {input.fq_trimmed_twice} | gzip > {output.fq_gz}
        """
# TODO, CHECK IF THE TRIMMING IS SUCCESSFUL, CHECK CROSS CONTAMINATION
rule fastqc_post_trim:
    input:
        "{libname}/fastqs/{sample_label}.umi.fqTrTr.rev.sorted.fq.gz"
    output:
        # PRPF8.umi.fqTrTr.rev.sorted_fastqc
        html="{libname}/fastqc/{sample_label}.umi.fqTrTr.rev.sorted_fastqc.html",
        txt="{libname}/fastqc/{sample_label}.umi.fqTrTr.rev.sorted_fastqc/fastqc_data.txt"
    threads: 2
    params:
        outdir="{libname}/fastqc/",
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


rule gather_fastqc_report from QC as qc_fastqc with:
    input:
        expand("{libname}/fastqc/{sample_label}.umi.fqTrTr.rev.sorted_fastqc/fastqc_data.txt", libname = libs, sample_label = rbps)
    output:
        basic='QC/{libname}/fastQC_basic_summary.csv',
        passfail='QC/{libname}/fastQC_passfail.csv'
    
use rule align_reads_to_Drosophila from snakeDros with:
    input:
        fq_1 = "{libname}/fastqs/{sample_label}.umi.fqTrTr.rev.sorted.fq.gz"
    output:
        ubam = "{libname}/bams/dros/{sample_label}.Aligned.out.bam",
        unmqpped= "{libname}/bams/dros/{sample_label}.Unmapped.out.mate1",
        log= "{libname}/bams/dros/{sample_label}.Log.final.out",
    params:
        error_out_file = "error_files/{sample_label}.{libname}_align_dros_reads",
        run_time = "08:40:00",
        cores = "4",
        memory = "10000",
        job_name = "align_reads",
        star_sjdb = STAR_DROS,
        outprefix = "{libname}/bams/dros/{sample_label}."
    benchmark: "benchmarks/align/{sample_label}.{libname}.align_dros_reads.txt"

use rule gather_mapstat from QC as mapstat_gather_repeat with:
    input:
        expand("{libname}/bams/dros/{sample_label}.Log.final.out", libname = libs, sample_label = rbps)
    output:
        "QC/dros_mapping_stats.csv"

rule align_reads_to_REPEAT:
    input:
        fq_1 = "{libname}/fastqs/{sample_label}.umi.fqTrTr.rev.sorted.fq.gz"
    output:
        ubam = "{libname}/bams/repeat/{sample_label}.Aligned.out.bam",
        unmqpped= "{libname}/bams/repeat/{sample_label}.Unmapped.out.mate1",
        log= "{libname}/bams/repeat/{sample_label}.Log.final.out",
    params:
        error_out_file = "error_files/{sample_label}_align_reads",
        run_time = "06:40:00",
        cores = "4",
        memory = "10000",
        job_name = "align_reads",
        star_sjdb = STAR_REP,
        outprefix = "{libname}/bams/repeat/{sample_label}.",
    benchmark: "benchmarks/align/{libname}.{sample_label}.align_reads.txt"
    shell:
        """
        module load star ;
        STAR \
            --alignEndsType EndToEnd \
            --genomeDir {params.star_sjdb} \
            --genomeLoad NoSharedMemory \
            --outBAMcompression 10 \
            --outFileNamePrefix {params.outprefix} \
            --outFilterMultimapNmax 30 \
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
use rule gather_mapstat from QC as mapstat_gather_repeat with:
    input:
        #find_all_files("{libname}/bams/repeat/{sample_label}.Log.final.out", libs)
        expand("{libname}/bams/repeat/{sample_label}.Log.final.out", libname = libs, sample_label = rbps)
    output:
        "QC/repeat_mapping_stats.csv"

rule align_to_GENOME:
    input:
        fq= "{libname}/bams/repeat/{sample_label}.Unmapped.out.mate1",
    output:
        ubam = "{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.out.bam",
        unmqpped= "{libname}/bams/genome/{sample_label}.genome-mapped.Unmapped.out.mate1",
        log= "{libname}/bams/genome/{sample_label}.genome-mapped.Log.final.out",
    params:
        error_out_file = "error_files/{libname}.{sample_label}_align_reads_genome",
        run_time = "06:40:00",
        cores = "4",
        memory = "10000",
        job_name = "align_reads",
        star_sjdb = STAR_DIR,
        outprefix = "{libname}/bams/genome/{sample_label}.genome-mapped.",
    benchmark: "benchmarks/align/{libname}.{sample_label}.align_reads.txt"
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

use rule gather_mapstat from QC as mapstat_gather_repeat with:
    input:
        #find_all_files("{libname}/bams/genome/{sample_label}.genome-mapped.Log.final.out", libs)
        expand("{libname}/bams/genome/{sample_label}.genome-mapped.Log.final.out", libname = libs, sample_label = rbps)
    output:
        "QC/genome_mapping_stats.csv"

rule sort_bams:
    input:
        bam="{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.out.bam",
    output:
        sort_once = "{libname}/bams/genome/{sample_label}.genome-mappedSo.Aligned.out.bam",
        sort_twice= "{libname}/bams/genome/{sample_label}.genome-mappedSoSo.Aligned.out.bam",
        bai = "{libname}/bams/genome/{sample_label}.genome-mappedSoSo.Aligned.out.bam.bai"
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
        bam="{libname}/bams/genome/{sample_label}.genome-mappedSoSo.Aligned.out.bam",
        bai="{libname}/bams/genome/{sample_label}.genome-mappedSoSo.Aligned.out.bam.bai"
    output:
        bam_dedup="{libname}/bams/genome/{sample_label}.genome-mappedSoSo.rmDup.Aligned.out.bam"
    params:
        error_out_file = "error_files/{libname}.{sample_label}_sort_bam",
        run_time = "06:40:00",
        cores = "4",
        memory = "10000",
        job_name = "sortbam",
        prefix='{libname}/bams/genome/{sample_label}.genome-mappedSoSo'
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
        bam = "{libname}/bams/genome/{sample_label}.genome-mappedSoSo.rmDup.Aligned.out.bam"
    output:
        sbam="{libname}/bams/genome/{sample_label}.genome-mappedSoSo.rmDupSo.Aligned.out.bam",
        bai = "{libname}/bams/genome/{sample_label}.genome-mappedSoSo.rmDupSo.Aligned.out.bam.bai"
    params:
        error_out_file = "error_files/{libname}.{sample_label}_index_bams",
        run_time = "01:40:00",
        cores = "4",
        memory = "1000",
        job_name = "index_bam"
    benchmark: "benchmarks/align/{libname}.{sample_label}.index_bam.txt"
    shell:
        "module load samtools;"
        "samtools sort -o {output.sbam} {input.bam} ;"
        "samtools index {output.sbam};"

use rule duplication_rate from QC as qc_duplication_rate with:
    input:
        dup=expand("{libname}/bams/genome/{sample_label}.genome-mappedSoSo.Aligned.out.bam", libname = libs, sample_label = rbps),
        rmdup=expand("{libname}/bams/genome/{sample_label}.genome-mappedSoSo.rmDupSo.Aligned.out.bam",libname = libs, sample_label = rbps)
    output:
        "QC/dup_level.csv"