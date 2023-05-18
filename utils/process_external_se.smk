# this file serves for users to process external control libraries the skipper way.
import pandas as pd
import glob
fastqs = ['/home/hsher/seqdata/20230401_charl_encode_totalrnaseq/ENCFF133KON.fastq.gz']
names = ['ENCODE_K562_totalRNAseq']
workdir:'/home/hsher/seqdata/20230401_charl_encode_totalrnaseq'
manifest = pd.DataFrame([fastqs, names], index = ['fastq1', 'Sample']).T
adapter_seq = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'

#snakemake -s utils/process_external_se.smk -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile /home/hsher/projects/oligoCLIP/config/preprocess_config/oligose_k562.yaml --use-conda --conda-prefix /home/hsher/snakeconda --use-singularity --singularity-prefix /home/hsher/singularity  -n
rule all:
    input:
        expand("output/bams/{sample_label}.Aligned.sortedByCoord.out.bam", sample_label = manifest['Sample'].tolist())

rule cutadapt:
    input:
        fq1=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["fastq1"].values[0]),
    output:
        fq1="output/fastqs/{sample_label}.Tr.fq.gz",
        metrics = "output/fastqs/{sample_label}.Tr.metrics",
    params:
        error_out_file = "error_files/cutadapt.{sample_label}",
        out_file = "stdout/cutadapt.{sample_label}",
        run_time = "6:45:00",
        cores = "4",
        memory = "10000",
        job_name = "cutadapt",
        adaptor1 = adapter_seq,
    benchmark: "benchmarks/cutadapt/{sample_label}.extract.txt"
    conda:
        "/home/hsher/projects/oligoCLIP/rules/envs/cutadapt.yaml"
    shell:
        """
        cutadapt \
            -a {params.adaptor1} \
            -o {output.fq1} \
            --times=2 \
            -m 23 \
            {input.fq1} > {output.metrics}
        """

rule align_reads:
    input:
        fq1="output/fastqs/{sample_label}.Tr.fq.gz",
    output:
        bam = "output/bams/{sample_label}.Aligned.sortedByCoord.out.bam",
        bai = "output/bams/{sample_label}.Aligned.sortedByCoord.out.bam.bai",
        unmapped1= "output/bams/{sample_label}.Unmapped.out.mate1",
        log= "output/bams/{sample_label}.Log.final.out",
    params:
        error_out_file = "error_files/{sample_label}.align_reads_genome.err",
        out_file = "stdout/{sample_label}.align_reads_genome.out",
        run_time = "02:00:00",
        memory = "40000",
        job_name = "align_reads",
        star_sjdb = config['STAR_DIR'],
        outprefix = "output/bams/{sample_label}.",
        cores = "8",
    benchmark: "benchmarks/align/{sample_label}.align_reads_genome.txt"
    shell:        
        """
        module load star
        STAR \
            --alignEndsType EndToEnd \
            --genomeDir {params.star_sjdb} \
            --genomeLoad NoSharedMemory \
            --outBAMcompression 10 \
            --outFileNamePrefix {params.outprefix} \
            --winAnchorMultimapNmax 100 \
            --outFilterMultimapNmax 100 \
            --outFilterMultimapScoreRange 1 \
            --outSAMmultNmax 1 \
            --outMultimapperOrder Random \
            --outFilterScoreMin 10 \
            --outFilterType BySJout \
            --limitOutSJcollapsed 5000000 \
            --outReadsUnmapped Fastx \
            --outSAMattrRGline ID:{wildcards.sample_label} \
            --outSAMattributes All \
            --outSAMmode Full \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outStd Log \
            --readFilesIn {input.fq1} \
            --readFilesCommand zcat \
            --runMode alignReads \
            --runThreadN {params.cores}\

        module load samtools
        samtools index {output.bam}
        """