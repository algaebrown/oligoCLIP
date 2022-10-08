QC_PATH = config['QC_TOOLS_PATH']

sample_labels = config['rbps']
libnames = config['libnames']
rbps = config['rbps']

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
        error_out_file = "error_files/qctrim.txt",
        out_file = "stdout/qctrim",
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
        error_out_file = "error_files/fastqc_stat.txt",
        out_file = "stdout/fastqc",
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
        out_file = "stdout/mapstat",
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
        out_file = "stdout/dup_rate",
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
        fq2=expand("{libname}/fastqs/ultraplex_demux_{sample_label}_Fwd.fastq.gz", libname = libnames, sample_label = rbps),
        fq1=expand("{libname}/fastqs/ultraplex_demux_{sample_label}_Rev.fastq.gz", libname = libnames, sample_label = rbps)
    output:
        'QC/demux_read_count.txt'
    params:
        error_out_file = "error_files/demux_count",
        out_file = "stdout/readcount",
        run_time = "00:40:00",
        cores = "1",
        memory = "10000",
        QC_PATH = QC_PATH,
    shell:
        """
        touch {output}
        for f in {input.fq2} ; do echo "$f $(zcat $f | wc -l)" >> {output}; done
        for f in {input.fq1} ; do echo "$f $(zcat $f | wc -l)" >> {output}; done
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
        out_file = "stdout/readwobar",
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

rule blast_unmapped_reads:
    input:
        target=HUMAN_RNA_NUCLEOTIDE,
        target_db=HUMAN_RNA_NUCLEOTIDE + '.nog',
        unmapped1_fq= "{libname}/bams/genome/{sample_label}.genome-mapped.Unmapped.out.mate1",
        unmapped2_fq= "{libname}/bams/genome/{sample_label}.genome-mapped.Unmapped.out.mate2"
    output:
        blast_result1="QC/unmapped_blast_output/{libname}.{sample_label}.1.blast.tsv",
        blast_result2="QC/unmapped_blast_output/{libname}.{sample_label}.2.blast.tsv",
        fasta1="QC/unmapped_blast_output/{libname}.{sample_label}.1.fasta",
        fasta2="QC/unmapped_blast_output/{libname}.{sample_label}.2.fasta"
    params:
        error_out_file = "error_files/demux_count",
        out_file = "stdout/blastunmap",
        run_time = "00:40:00",
        cores = "1",
        nlines = N_READ_TO_SAMPLE * 4
    shell:
        """
        set +o pipefail; 
        module load samtools
        samtools fasta {input.unmapped1_fq} | head -n {params.nlines} > {output.fasta1}
        samtools fasta {input.unmapped2_fq} | head -n {params.nlines} > {output.fasta2}
        module load blast
        blastn -db {input.target} -query {output.fasta1} -out {output.blast_result1} -outfmt 6 -max_target_seqs 1 
        blastn -db {input.target} -query {output.fasta2} -out {output.blast_result2} -outfmt 6 -max_target_seqs 1 
        """

rule blast_unmapped_reads_too_short:
    input:
        target=HUMAN_RNA_NUCLEOTIDE,
        target_db=HUMAN_RNA_NUCLEOTIDE + '.nog',
        bam= "{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam",
    output:
        blast_result="QC/unmapped_blast_output/{libname}.{sample_label}.short.blast.tsv",
        fasta="QC/unmapped_blast_output/{libname}.{sample_label}.short.fasta",
    params:
        error_out_file = "error_files/demux_count",
        out_file = "stdout/blastunmap",
        run_time = "00:40:00",
        cores = "1",
        nlines = N_READ_TO_SAMPLE
    shell:
        """
        set +o pipefail; 
        module load samtools
        samtools view -f 4 {input.bam} | grep uT:A:1 | head -n {params.nlines} | samtools fasta >  {output.fasta}
        module load blast
        blastn -db {input.target} -query {output.fasta} -out {output.blast_result} -outfmt 6 -max_target_seqs 1 
        """