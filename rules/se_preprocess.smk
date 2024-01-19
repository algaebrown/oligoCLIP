import pandas as pd
locals().update(config)

manifest = pd.read_csv(config['MANIFEST'])
barcode_df = pd.read_csv(config['barcode_csv'], header = None, sep = ':', names = ['barcode', 'RBP'])
rbps = barcode_df['RBP'].tolist()

rule fastqc_initial:
    input:
        fq1 = lambda wildcards: manifest.loc[manifest['libname']==wildcards.libname, 'fastq'],
    output:
        html1="{libname}/fastqc/initial_r1_fastqc.html",
        txt1="{libname}/fastqc/initial_r1_fastqc/fastqc_data.txt",
    params:
        outdir="{libname}/fastqc/",
        run_time = "02:09:00",
        cores="2",
        error_out_file = "error_files/fastqc.{libname}.txt",
        out_file = "stdout/fastqc.{libname}.txt",
        memory = 20000,
    benchmark: "benchmarks/qc/fastqc_initial.{libname}.txt"
    container:
        "docker://howardxu520/skipper:fastqc_0.12.1"
    shell:
        """
        zcat {input.fq1} | \
            fastqc stdin:initial_r1 \
            --extract \
            --outdir {params.outdir} \
            -t {params.cores}
        """

rule tile_adaptor:
    output:
        adafwd = temp("params/adaptor_fwd.fasta"),
    params:
        run_time = "00:20:00",
        cores="1",
        adaptor_fwd = config['adaptor_fwd'],
        error_out_file = "error_files/adatile.txt",
        out_file = "stdout/adatile.txt",
        tiling_length = config['tile_length'],
        memory = 10000,
    conda:
        "envs/metadensity.yaml"
    benchmark: "benchmarks/tile_adaptor"
    shell:      
        """
        python {SCRIPT_PATH}/create_adaptor_tile.py {params.adaptor_fwd} {output.adafwd} {params.tiling_length}
        """

rule trim_adaptor:
    input:
        fq = lambda wildcards: manifest.loc[manifest['libname']==wildcards.libname, 'fastq'].iloc[0],
        ada_fwd = rules.tile_adaptor.output.adafwd
    output:
        fq_trimmed=temp("{libname}/fastqs/all.Tr.fq.gz"),
        metrics = "QC/{libname}.Tr.metrics"
    params:
        run_time = "12:04:00",
        cores="4",
        error_out_file = "error_files/trim_adaptor.{libname}.txt",
        out_file = "stdout/trim_adaptor.{libname}.txt",
        quality_cutoff = config['QUALITY_CUTOFF'],
        memory = 80000,
    benchmark: "benchmarks/cutadapt/trim_adaptor.{libname}.txt"
    conda:
        "envs/cutadapt.yaml"
    shell:
        """
        cutadapt -O 1 \
            --times 2 \
            -e 0.1 \
            --quality-cutoff {params.quality_cutoff} \
            -m 23 \
            -o {output.fq_trimmed} \
            -a file:{input.ada_fwd} \
            --cores=0 \
            {input.fq} > {output.metrics}
        """

rule extract_umi: 
    input:
        fq = rules.trim_adaptor.output.fq_trimmed
    output:
        fq_umi = temp("{libname}/fastqs/all.Tr.umi.fq.gz"),
        metrics = "QC/{libname}.all.Tr.umi.metrics"
    params:
        error_out_file = "error_files/extract_umi",
        out_file = "stdout/extract_umi.{libname}.txt",
        run_time = "3:45:00",
        cores = "4",
        umi_pattern = config['umi_pattern'],
        memory = 80000,
    benchmark: "benchmarks/umi/extract.{libname}.txt"
    conda:
        "envs/umi_tools.yaml"
    shell:
        """
        umi_tools extract \
            --random-seed 1 \
            --bc-pattern {params.umi_pattern} \
            --stdin {input.fq} \
            --stdout {output.fq_umi} \
            --log {output.metrics} \
        """

rule demultiplex:
    input:
        fq1 = rules.extract_umi.output.fq_umi,
        barcode_csv = ancient(config['barcode_csv'])
    output:
        fq=temp(expand("{libname}/fastqs/ultraplex_demux_{sample_label}.fastq.gz", 
            libname = ["{libname}"], sample_label = rbps)),
        missing_fq="{libname}/fastqs/ultraplex_demux_5bc_no_match.fastq.gz",
        # logs = "{libname}/barcode.log"
    conda:
        "envs/ultraplex.yaml"
    benchmark: "benchmarks/pre/demux.{libname}.txt"
    params:
        outdir="",
        cores="4",
        run_time = "06:00:00",
        prefix = "{libname}/fastqs/",
        error_out_file = "error_files/demux.{libname}.txt",
        out_file = "stdout/demux.{libname}.txt",
        memory = 160000,
    shell:
        """
        cd {params.prefix}
        # in case there is no not matched ones
        
        ultraplex -i all.Tr.umi.fq.gz -b {input.barcode_csv}  \
            -m5 1 -m3 0 -t {params.cores} -a XX -a2 XX --ultra
        if [ ! -f ultraplex_demux_5bc_no_match.fastq.gz ]
        then
            touch ultraplex_demux_5bc_no_match.fastq.gz
        fi
        """

rule reverse_complement:
    input:
        "{libname}/fastqs/ultraplex_demux_{sample_label}.fastq.gz"
    output:
        fqrev="{libname}/fastqs/ultraplex_demux_{sample_label}.rev.fastq.gz"
    params:
        run_time = "05:04:00",
        cores="1",
        out_file = "stdout/revcomp.{libname}.txt",
        error_out_file = "error_files/revcomp.{libname}.txt",
        memory = 40000,
    conda:
        "envs/fastx_toolkit.yaml"
    shell:
        """
        zcat {input} | fastx_reverse_complement | gzip > {output.fqrev}
        """

rule fastqc_post_trim:
    input:
        fq1=rules.reverse_complement.output.fqrev
    output:
        html1="{libname}/fastqc/ultraplex_demux_{sample_label}.rev_fastqc.html",
        txt1="{libname}/fastqc/ultraplex_demux_{sample_label}.rev_fastqc/fastqc_data.txt",
    params:
        outdir="{libname}/fastqc/",
        run_time = "02:09:00",
        cores="1",
        error_out_file = "error_files/fastqc.{libname}.txt",
        out_file = "stdout/fastqc.{libname}.txt",
        memory = 40000,
    benchmark: "benchmarks/qc/fastqc.{libname}.{sample_label}.txt"
    container:
        "docker://howardxu520/skipper:fastqc_0.12.1"
    shell:
        """
        fastqc {input.fq1} --extract --outdir {params.outdir} -t {params.cores}
        """

rule align_reads:
    input:
        fq1=rules.reverse_complement.output.fqrev
    output:
        bam = temp("{libname}/bams/{sample_label}.Aligned.sortedByCoord.out.bam"),
        unmapped1= "{libname}/bams/{sample_label}.Unmapped.out.mate1",
        log= "{libname}/bams/{sample_label}.Log.final.out",
    params:
        error_out_file = "error_files/{libname}.{sample_label}.align_reads_genome.err",
        out_file = "stdout/{libname}.{sample_label}.align_reads_genome.out",
        run_time = "02:00:00",
        star_sjdb = config['STAR_DIR'],
        outprefix = "{libname}/bams/{sample_label}.",
        cores = "8",
        memory = 320000,
    benchmark: "benchmarks/align/{libname}.{sample_label}.align_reads_genome.txt"
    container:
        "docker://howardxu520/skipper:star_2.7.10b"
    shell:        
        """
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
            --runThreadN {params.cores}
        """

rule umi_dedup:
    input:
        bam=rules.align_reads.output.bam,
        bai="{libname}/bams/{sample_label}.Aligned.sortedByCoord.out.bam.bai"
    output:
        bam_dedup="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
    params:
        error_out_file = "error_files/dedup.{libname}.{sample_label}",
        out_file = "stdout/{libname}.{sample_label}.index_reads",
        run_time = "06:40:00",
        cores = "4",
        prefix='{libname}/bams/genome/{sample_label}.genome-mapped',
        memory = 160000,
    benchmark: "benchmarks/align/dedup.{libname}.{sample_label}.txt"
    container:
        "docker://howardxu520/skipper:umicollapse_1.0.0"
    shell:
        """
        java -server -Xms8G -Xmx8G -Xss20M -jar /UMICollapse/umicollapse.jar bam -i {input.bam} -o {output.bam_dedup} --umi-sep : --two-pass --paired
        """

rule index_bam:
    input:
        "{anything}.bam"
    output:
        "{anything}.bam.bai"
    params:
        error_out_file = "error_files/index_bam",
        out_file = "stdout/index_bam",
        run_time = "40:00",
        cores = "1",
        memory = 40000,
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools index {input}
        """