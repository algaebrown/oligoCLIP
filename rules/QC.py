QC_PATH = config['QC_TOOLS_PATH']

sample_labels = []
libnames = []
rbps = []

HUMAN_RNA_NUCLEOTIDE='/projects/ps-yeolab4/seqdata/20200622_gencode_coords_hsher/GRCh38_latest_rna.fna'
N_READ_TO_SAMPLE=5*10**3

rule gather_trimming_stat:
    input:
        tr1=expand("fastqs/umi/{sample_label}.umi.r1.fqTr.metrics", sample_label = sample_labels)
    output:
        tr1='QC/cutadapt_log1.csv',
    params:
        run_time = "00:10:00",
        cores="1",
        QC_PATH = QC_PATH,
        error_out_file = "error_files/qctrim.txt"
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {params.QC_PATH}/trimming_stat.py "{input.tr1:q}" {output.tr1}
        """

rule gather_fastqc_report:
    input:
        expand("fastqc/{sample_label}.umi.fqTrTr_fastqc/fastqc_data.txt", sample_label = sample_labels)
    output:
        basic='QC/fastQC_basic_summary.csv',
        passfail='QC/fastQC_passfail.csv'
    params:
        run_time = "00:40:00",
        cores = "1",
        memory = "10000",
        job_name = "gather_stat",
        QC_PATH = QC_PATH,
        error_out_file = "error_files/fastqc_stat.txt"
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {params.QC_PATH}/fastqc_io.py -i "{input}" -p {output.passfail} -b {output.basic}
        """

rule gather_mapstat:
    input:
        #find_all_files("{libname}/bams/repeat/{sample_label}.Log.final.out", libs)
        expand("bams/dros/{sample_label}.Log.final.out", sample_label = sample_labels)
    output:
        "QC/dros_mapping_stats.csv"
    conda:
        "envs/metadensity.yaml"
    params:
        error_out_file = "error_files/mapstat",
        run_time = "00:40:00",
        cores = "1",
        memory = "10000",
        job_name = "gather_stat",
        QC_PATH = QC_PATH,
    shell:
        """
        python {params.QC_PATH}/star_mapping_stat_io.py -i "{input}" -o {output}
        """

rule duplication_rate:
    input:
        dup=expand("{libname}/bams/genome/{sample_label}.genome-mapped.Log.final.out", libname = libnames, sample_label = rbps),
        rmdup=expand("{libname}/bams/genome/{sample_label}.genome-mapped.rmDup.Aligned.sortedByCoord.out.bam",libname = libnames, sample_label = rbps)
    output:
        'QC/dup_level.csv'
    params:
        error_out_file = "error_files/dup",
        run_time = "00:40:00",
        cores = "1",
        memory = "10000",
        job_name = "gather_stat",
        QC_PATH = QC_PATH,
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {params.QC_PATH}/dup_level.py "{input.dup}" "{input.rmdup}" {output}
        """

rule count_demultiplex_ultraplex:
    input:
        fq2=expand("{libname}/fastqs/ultraplex_demux_{sample_label}_Fwd.fastq.gz", libname = libnames, sample_label = rbps)
    output:
        'QC/demux_read_count.txt'
    params:
        error_out_file = "error_files/demux_count",
        run_time = "00:40:00",
        cores = "1",
        memory = "10000",
        QC_PATH = QC_PATH,
    shell:
        """
        wc -l {input} > {output}
        """

rule what_is_read_wo_barcode:
    input:
        target=HUMAN_RNA_NUCLEOTIDE,
        target_db=HUMAN_RNA_NUCLEOTIDE + '.nog',
        query_fq_gz="{libname}/fastqs/ultraplex_demux_5bc_no_match_Fwd.fastq.gz"
    output:
        blast_result="QC/nobarcode_blast_output/{libname}.blast.tsv",
        fasta="QC/nobarcode_blast_output/{libname}.fasta"
    params:
        error_out_file = "error_files/demux_count",
        run_time = "00:40:00",
        cores = "1",
        nlines = N_READ_TO_SAMPLE * 4
    shell:
        """
        set +o pipefail; 
        zcat {input.query_fq_gz} | head -n {params.nlines} | sed -n '1~4s/^@/>/p;2~4p' > {output.fasta}
        module load blast
        blastn -db {input.target} -query {output.fasta} -out {output.blast_result} -outfmt 6 -max_target_seqs 1 
        """