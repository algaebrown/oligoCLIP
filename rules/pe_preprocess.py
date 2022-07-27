import pandas as pd
manifest = pd.read_csv(config['MANIFEST'])
# /projects/ps-yeolab3/eboyle/encode/pipeline/04_20220506

barcode_df = pd.read_csv(config['barcode_csv'], header = None, sep = ':', names = ['barcode', 'RBP'])
rbps = barcode_df['RBP'].tolist()
rule trim_adaptor:
    input:
        fq1 = lambda wildcards: manifest.loc[manifest['libname']==wildcards.libname, 'fastq1'],
        fq2 = lambda wildcards: manifest.loc[manifest['libname']==wildcards.libname, 'fastq2']
    output:
        fq1 = "{libname}/fastqs/all.Tr.fq1.gz",
        fq2 = "{libname}/fastqs/all.Tr.fq2.gz",
        metric = "QC/{libname}.Tr.metrics"
    params:
        run_time = "12:04:00",
        cores="4",
        adaptor_fwd = config['adaptor_fwd'],
        adaptor_rev = config['adaptor_rev'],
        error_out_file = "error_files/trim_adaptor.{libname}.txt",
        quality_cutoff = config['QUALITY_CUTOFF']
    conda:
        "envs/cutadapt.yaml"
    benchmark: "benchmarks/trim_adaptor.{libname}"
    shell:      
        """
        cutadapt -a {params.adaptor_fwd} \
            -A {params.adaptor_rev} \
            --times 2 \
            -e 0.1 \
            --quality-cutoff {params.quality_cutoff} \
            -m 23 \
            -o {output.fq1} \
            -p {output.fq2} \
            --cores=0 \
            {input.fq1} {input.fq2} > {output.metric}
        """

rule extract_umi_and_trim_polyG: # TODO: adaptor TRIM first
    input:
        fq1 = "{libname}/fastqs/all.Tr.fq1.gz",
        fq2 = "{libname}/fastqs/all.Tr.fq2.gz",
    output:
        fq1 = "{libname}/fastqs/all.Tr.umi.fq1.gz",
        fq2 = "{libname}/fastqs/all.Tr.umi.fq2.gz",
        metrics = "QC/{libname}.umi.json",
        metrics2 = "QC/{libname}.umi.html",
    params:
        error_out_file = "error_files/extract_umi.{libname}.txt",
        run_time = "3:45:00",
        cores = "4",
        memory = "10000",
        job_name = "extract_umi",
        umi_length = config['umi_length']
    conda:
        "envs/fastp.yaml"
    benchmark: "benchmarks/umi/extract_umi.{libname}.txt"
    shell:
        """
        fastp -i {input.fq1} -I {input.fq2} \
            -o {output.fq1} -O {output.fq2} \
            --disable_adapter_trimming \
            --umi \
            --umi_len={params.umi_length} \
            --umi_loc=read1 \
            --trim_poly_g \
            -j {output.metrics} \
            -h {output.metrics2} \
            -w {params.cores}
        """

# reverse read1 and read2 cause ultraplex does not support 3' only demux
# set adaptor to X to disable adaptor trimming
# the pgzip thing is slow. always lead to incomplete file
rule demultiplex:
    input:
        fq1 = "{libname}/fastqs/all.Tr.umi.fq1.gz",
        fq2 = "{libname}/fastqs/all.Tr.umi.fq2.gz",
        barcode_csv = ancient(config['barcode_csv'])
    output:
        fq1=expand("{libname}/fastqs/ultraplex_demux_{sample_label}_Rev.fastq.gz", libname = ["{libname}"], sample_label = rbps),
        fq2=expand("{libname}/fastqs/ultraplex_demux_{sample_label}_Fwd.fastq.gz", libname = ["{libname}"], sample_label = rbps),
        missing_fq1="{libname}/fastqs/ultraplex_demux_5bc_no_match_Rev.fastq.gz",
        missing_fq2="{libname}/fastqs/ultraplex_demux_5bc_no_match_Fwd.fastq.gz",
        # logs = "{libname}/barcode.log"
    conda:
        "envs/ultraplex.yaml"
    benchmark: "benchmarks/pre/demux.{libname}.txt"
    params:
        outdir="",
        cores="4",
        run_time = "06:00:00",
        prefix = "{libname}/fastqs/",
        error_out_file = "error_files/demux.{libname}.txt"
    shell:
        """
        cd {params.prefix}
        ultraplex -i all.Tr.umi.fq2.gz -i2 all.Tr.umi.fq1.gz -b {input.barcode_csv}  \
            -m5 1 -m3 0 -t {params.cores} -a XX -a2 XX --ultra
        sleep 600
        """

# TODO, CHECK IF THE TRIMMING IS SUCCESSFUL, CHECK CROSS CONTAMINATION
rule fastqc_post_trim:
    input:
        fq1="{libname}/fastqs/ultraplex_demux_{sample_label}_Rev.fastq.gz",
        fq2="{libname}/fastqs/ultraplex_demux_{sample_label}_Fwd.fastq.gz"
    output:
        # PRPF8.umi.fqTrTr.rev.sorted_fastqc
        html1="{libname}/fastqc/ultraplex_demux_{sample_label}_Rev_fastqc.html",
        txt1="{libname}/fastqc/ultraplex_demux_{sample_label}_Rev_fastqc/fastqc_data.txt",
        html2="{libname}/fastqc/ultraplex_demux_{sample_label}_Fwd_fastqc.html",
        txt2="{libname}/fastqc/ultraplex_demux_{sample_label}_Fwd_fastqc/fastqc_data.txt"
    params:
        outdir="{libname}/fastqc/",
        run_time = "02:09:00",
        cores="1",
        error_out_file = "error_files/fastqc.{libname}.txt"
    benchmark: "benchmarks/qc/fastqc.{libname}.{sample_label}.txt"
    shell:
        """
        module load fastqc;
        fastqc {input.fq1} --extract --outdir {params.outdir} -t {params.cores}
        fastqc {input.fq2} --extract --outdir {params.outdir} -t {params.cores}
        """


rule align_reads_to_REPEAT:
    input:
        fq1="{libname}/fastqs/ultraplex_demux_{sample_label}_Rev.fastq.gz",
        fq2="{libname}/fastqs/ultraplex_demux_{sample_label}_Fwd.fastq.gz"
    output:
        ubam = "{libname}/bams/repeat/{sample_label}.Aligned.out.bam",
        unmapped1= "{libname}/bams/repeat/{sample_label}.Unmapped.out.mate1",
        unmapped2= "{libname}/bams/repeat/{sample_label}.Unmapped.out.mate2",
        log= "{libname}/bams/repeat/{sample_label}.Log.final.out",
    params:
        error_out_file = "error_files/{sample_label}_align_reads",
        run_time = "06:40:00",
        cores = "4",
        memory = "10000",
        job_name = "align_reads",
        star_sjdb = config['STAR_REP'],
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
            --readFilesIn {input.fq1} {input.fq2}\
            --readFilesCommand zcat \
            --runMode alignReads \
            --runThreadN 8
        """

rule align_to_GENOME:
    input:
        fq1= "{libname}/bams/repeat/{sample_label}.Unmapped.out.mate1",
        fq2= "{libname}/bams/repeat/{sample_label}.Unmapped.out.mate2",
    output:
        bam = "{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam",
        bai = "{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam.bai",
        unmapped1= "{libname}/bams/genome/{sample_label}.genome-mapped.Unmapped.out.mate1",
        unmapped2= "{libname}/bams/genome/{sample_label}.genome-mapped.Unmapped.out.mate2",
        log= "{libname}/bams/genome/{sample_label}.genome-mapped.Log.final.out",
    params:
        error_out_file = "error_files/{libname}.{sample_label}_align_reads_genome",
        run_time = "06:40:00",
        cores = "4",
        memory = "10000",
        job_name = "align_reads",
        star_sjdb = config['STAR_DIR'],
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
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outStd Log \
        --readFilesIn {input.fq1} {input.fq2} \
        --runMode alignReads \
        --runThreadN 8

        module load samtools
        samtools index {output.bam}
        """


rule umi_dedup:
    # TODO: containerize
    input:
        bam="{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam",
        bai="{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam.bai",
    output:
        bam_dedup="{libname}/bams/genome/{sample_label}.genome-mapped.rmDup.Aligned.sortedByCoord.out.bam"
    params:
        error_out_file = "error_files/dedup.{libname}.{sample_label}",
        run_time = "06:40:00",
        cores = "4",
        memory = "10000",
        job_name = "sortbam",
        prefix='{libname}/bams/genome/{sample_label}.genome-mapped',
        java = config['JAVA_PATH'],
        umicollapse = config['UMICOLLAPSE_PATH']
    conda: 
        "envs/umi_tools.yaml"
    benchmark: "benchmarks/align/dedup.{libname}.{sample_label}.txt"
    shell:
        """
        {params.java} -server -Xms8G -Xmx8G -Xss20M -jar {params.umicollapse} bam -i {input.bam} -o {output.bam_dedup} --umi-sep : --two-pass --paired
        """
rule index_genome_bams:
    input:
        bam = "{libname}/bams/genome/{sample_label}.genome-mapped.rmDup.Aligned.sortedByCoord.out.bam"
    output:
        bai = "{libname}/bams/genome/{sample_label}.genome-mapped.rmDup.Aligned.sortedByCoord.out.bam.bai"
    params:
        error_out_file = "error_files/index_bam.{libname}.{sample_label}",
        run_time = "01:40:00",
        cores = "4",
        memory = "1000",
        job_name = "index_bam"
    benchmark: "benchmarks/align/{libname}.{sample_label}.index_bam.txt"
    shell:
        "module load samtools;"
        "samtools index {input.bam};"
