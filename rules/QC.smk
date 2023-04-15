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
        dup=expand("{libname}/bams/{sample_label}.Log.final.out", libname = libnames, sample_label = rbps),
        rmdup=expand("{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",libname = libnames, sample_label = rbps),
        rmdup_bai=expand("{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam.bai",libname = libnames, sample_label = rbps)
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

#echo haha $(echo $(zcat multiplex_HEK293_3/fastqs/ultraplex_demux_QKI_Rev.fastq.gz | wc -l)/4 | bc)
rule count_demultiplex_ultraplex:
    input:
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
        for f in {input.fq1} ; do \
        echo "$f $(echo $(zcat $f | wc -l)/4 | bc)" >> {output}; \
        done
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
        unmapped1_fq= "{libname}/bams/{sample_label}.Unmapped.out.mate1",
        unmapped2_fq= "{libname}/bams/{sample_label}.Unmapped.out.mate2"
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
        bam= "{libname}/bams/{sample_label}.Aligned.sortedByCoord.out.bam",
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
rule make_read_count_summary:
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        counts = "counts/genome/megatables/{libname}.tsv.gz",
    output:
        region_summary =  "QC/read_count/{libname}.region.csv",
        type_summary =  "QC/read_count/{libname}.genetype.csv",
        name_summary =  "QC/read_count/{libname}.genename.csv",
        dist = "QC/read_count/{libname}.cosine_similarity.csv"
    params:
        error_out_file = "error_files/{libname}.read_count_summary.err",
        out_file = "stdout/{libname}.read_count_summary.out",
        run_time = "1:20:00",
        cores = 1
    run:
        import os
        print(output.region_summary)
        try:
            os.mkdir('QC/read_count')
        except Exception as e:
            print(e)
        import pandas as pd
        cnt = pd.read_csv(input.counts, sep = '\t')
        feature_annotations = pd.read_csv(input.feature_annotations, sep = '\t')

        df = pd.concat([feature_annotations, cnt], axis = 1)

        by_type = df.groupby(by = 'feature_type_top')[cnt.columns].sum()
        by_gene = df.groupby(by = 'gene_type_top')[cnt.columns].sum()
        by_name = df.groupby(by = 'gene_name')[cnt.columns].sum()

        by_type.to_csv(output.region_summary)
        by_gene.to_csv(output.type_summary)
        by_name.to_csv(output.name_summary)

        # distance
        from scipy.spatial.distance import pdist, squareform
        cov_filter = 10
        dist = squareform(pdist(cnt.loc[cnt.sum(axis = 1)>cov_filter].T, 'cosine'))

        dist_df = pd.DataFrame(1-dist, columns = cnt.columns, index = cnt.columns)
        dist_df.to_csv(output.dist)
