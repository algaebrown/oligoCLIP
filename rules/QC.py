rule gather_trimming_stat:
    input:
        tr1=expand("fastqs/umi/{sample_label}.umi.r1.fqTr.metrics", sample_label = sample_labels)
    output:
        tr1='QC/cutadapt_log1.csv',
    params:
        run_time = "00:10:00",
        cores="1",
        files2=','.join(expand("fastqs/umi/{sample_label}.umi.r1.fqTrTr.metrics", sample_label = sample_labels)),
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python /home/hsher/projects/QC_tools/trimming_stat.py {params.files1} {output.tr1}
        """

rule gather_fastqc_report:
    input:
        expand("fastqc/{sample_label}.umi.fqTrTr_fastqc/fastqc_data.txt", sample_label = sample_labels)
    output:
        basic='QC/fastQC_basic_summary.csv',
        passfail='QC/fastQC_passfail.csv'
    params:
        error_out_file = "error_files/mapstat",
        run_time = "00:40:00",
        cores = "1",
        memory = "10000",
        job_name = "gather_stat",
        files = ','.join(expand("fastqc/{sample_label}.umi.fqTrTr_fastqc/fastqc_data.txt", 
            sample_label = sample_labels))
    shell:
        """
        python /home/hsher/projects/QC_tools/fastqc_io.py -i {params.files} -p {output.passfail} -b {output.basic}
        """

rule gather_DROSOPHILA_mapstat:
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
        files = ','.join(expand("bams/dros/{sample_label}.Log.final.out", sample_label = sample_labels))
    shell:
        """
        python /home/hsher/projects/QC_tools/star_mapping_stat_io.py -i {params.files} -o {output}
        """