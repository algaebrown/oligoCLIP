import pandas as pd
manifest = pd.read_csv(config['MANIFEST'])
SCRIPT_PATH=config['SCRIPT_PATH']
# /projects/ps-yeolab3/eboyle/encode/pipeline/04_20220506

barcode_df = pd.read_csv(config['barcode_csv'], header = None, sep = ':', names = ['barcode', 'RBP'])
rbps = barcode_df['RBP'].tolist()

rule tile_adaptor:
    output:
        adafwd = "params/adaptor_fwd.fasta",
        adarev = "params/adaptor_rev.fasta"
    params:
        run_time = "00:20:00",
        cores="1",
        adaptor_fwd = config['adaptor_fwd'],
        adaptor_rev = config['adaptor_rev'],
        error_out_file = "error_files/adatile.txt",
        out_file = "stdout/adatile.txt",
        tiling_length = config['tile_length']
    conda:
        "envs/metadensity.yaml"
    benchmark: "benchmarks/tile_adaptor"
    shell:      
        """
        python {SCRIPT_PATH}create_adaptor_tile.py {params.adaptor_fwd} {output.adafwd} {params.tiling_length}
        python {SCRIPT_PATH}create_adaptor_tile.py {params.adaptor_rev} {output.adarev} {params.tiling_length}
        """

rule trim_adaptor:
    input:
        fq1 = lambda wildcards: manifest.loc[manifest['libname']==wildcards.libname, 'fastq1'],
        fq2 = lambda wildcards: manifest.loc[manifest['libname']==wildcards.libname, 'fastq2'],
        adaptor_fwd = "params/adaptor_fwd.fasta",
        adaptor_rev = "params/adaptor_rev.fasta"
    output:
        fq1 = "{libname}/fastqs/all.Tr.fq1.gz",
        fq2 = "{libname}/fastqs/all.Tr.fq2.gz",
        metric = "QC/{libname}.Tr.metrics"
    params:
        run_time = "12:04:00",
        cores="4",
        error_out_file = "error_files/trim_adaptor.{libname}.txt",
        quality_cutoff = config['QUALITY_CUTOFF'],
        out_file = "stdout/trim_adaptor.{libname}.txt.txt",
    conda:
        "envs/cutadapt.yaml"
    benchmark: "benchmarks/pre/trim_adaptor.{libname}"
    shell:      
        """
        cutadapt -a file:{input.adaptor_fwd} \
            -A file:{input.adaptor_rev} \
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
#fastp does not remove UMI from read2
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
        out_file = "stdout/extract_umi.{libname}.txt",
        run_time = "3:45:00",
        cores = "4",
        memory = "10000",
        job_name = "extract_umi",
        umi_length = config['umi_length']
    conda:
        "envs/fastp.yaml"
    benchmark: "benchmarks/pre/extract_umi.{libname}.txt"
    shell:
        """
        fastp -i {input.fq1} -I {input.fq2} \
            -o {output.fq1} -O {output.fq2} \
            --disable_adapter_trimming \
            --disable_length_filtering \
            --umi \
            --umi_len={params.umi_length} \
            --umi_loc=read1 \
            --trim_poly_g \
            -j {output.metrics} \
            -h {output.metrics2} \
            -w {params.cores}
        """

rule trim_umi_from_read2:
    # sometimes read2 has some UMI left despite fastp claims to trim it.
    input:
        fq2 = "{libname}/fastqs/all.Tr.umi.fq2.gz",
    output:
        fq2_unzip = temp("{libname}/fastqs/all.Tr.umi.fq2"),
        fq2 = "{libname}/fastqs/all.Tr.umi.fq2.trim.gz",
    params:
        error_out_file = "error_files/remove_UMI_from_r2.{libname}.txt",
        out_file = "stdout/remove_UMI_from_r2.{libname}.txt",
        run_time = "3:45:00",
        cores = "4",
        memory = "10000",
        job_name = "extract_umi",
        umi_length = config['umi_length']
    conda:
        "envs/seqtk.yaml"
    benchmark: "benchmarks/pre/trim_umi_from_r2.{libname}.txt"
    shell:
        """
        zcat {input.fq2} > {output.fq2_unzip}
        seqtk trimfq -e {params.umi_length} {output.fq2_unzip} | gzip > {output.fq2}
        """
# reverse read1 and read2 cause ultraplex does not support 3' only demux
# set adaptor to X to disable adaptor trimming
# the pgzip thing is slow. always lead to incomplete file
rule demultiplex: ################ never used the trimmed fastq file
    input:
        fq1 = "{libname}/fastqs/all.Tr.umi.fq1.gz",
        fq2 = "{libname}/fastqs/all.Tr.umi.fq2.trim.gz",
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
        error_out_file = "error_files/demux.{libname}.txt",
        out_file = "stdout/demux.{libname}.txt",
    shell:
        """
        cd {params.prefix}
        # in case there is no not matched ones
        
        ultraplex -i all.Tr.umi.fq2.gz -i2 all.Tr.umi.fq2.trim.gz -b {input.barcode_csv}  \
            -m5 1 -m3 0 -t {params.cores} -a XX -a2 XX --ultra
        if [ ! -f ultraplex_demux_5bc_no_match_Rev.fastq.gz ]
        then
            touch ultraplex_demux_5bc_no_match_Rev.fastq.gz
            touch ultraplex_demux_5bc_no_match_Fwd.fastq.gz
        fi
        """

rule trim_barcode_r1:
# ultraplex does not trim well for the other read
    input:
        fq1 = "{libname}/fastqs/ultraplex_demux_{sample_label}_Rev.fastq.gz", # has reverse complement of barcode at the end
        fq2 = "{libname}/fastqs/ultraplex_demux_{sample_label}_Fwd.fastq.gz",
    output:
        fq1 = "{libname}/fastqs/ultraplex_demux_{sample_label}_Rev.Tr.fastq.gz",
        fq2 = "{libname}/fastqs/ultraplex_demux_{sample_label}_Fwd.Tr.fastq.gz",
        metric = "QC/{libname}.{sample_label}.n_rev_read_w_bar.txt",
        metric_html = "QC/{libname}.{sample_label}.n_rev_read_w_bar.html",
        metric_json = "QC/{libname}.{sample_label}.n_rev_read_w_bar.json"
    params:
        run_time = "12:04:00",
        cores="4",
        error_out_file = "error_files/trim_barcode.{libname}.{sample_label}.txt",
        quality_cutoff = config['QUALITY_CUTOFF'],
        out_file = "stdout/trim_barcode.{libname}.{sample_label}.txt",
        barcode_sequence = lambda wildcards:barcode_df.loc[barcode_df['RBP'] == wildcards.sample_label,'barcode'].iloc[0],
    conda:
        "envs/fastp.yaml"
    benchmark: "benchmarks/pre/trim_barcode.{libname}.{sample_label}"
    shell:
        """
        set +o pipefail; 
        rev_bar=$(echo {params.barcode_sequence} | tr ACGTacgt TGCAtgca | rev)
        fastp -i {input.fq1} \
            -I {input.fq2} \
            -o {output.fq1} \
            -O {output.fq2} \
            --adapter_sequence $rev_bar \
            -j {output.metric_json} \
            -h {output.metric_html} \
            --disable_length_filtering
        zcat {input.fq2} | grep -v "@" | grep $rev_bar | wc -l > {output.metric}
        zcat {output.fq2} | grep -v "@" | grep $rev_bar | wc -l >> {output.metric}
        """


# TODO, CHECK IF THE TRIMMING IS SUCCESSFUL, CHECK CROSS CONTAMINATION
rule fastqc_post_trim:
    input:
        fq1="{libname}/fastqs/ultraplex_demux_{sample_label}_Rev.Tr.fastq.gz",
        fq2="{libname}/fastqs/ultraplex_demux_{sample_label}_Fwd.Tr.fastq.gz",
    output:
        # PRPF8.umi.fqTrTr.rev.sorted_fastqc
        html1="{libname}/fastqc/ultraplex_demux_{sample_label}_Rev.Tr_fastqc.html",
        txt1="{libname}/fastqc/ultraplex_demux_{sample_label}_Rev.Tr_fastqc/fastqc_data.txt",
        html2="{libname}/fastqc/ultraplex_demux_{sample_label}_Fwd.Tr_fastqc.html",
        txt2="{libname}/fastqc/ultraplex_demux_{sample_label}_Fwd.Tr_fastqc/fastqc_data.txt"
    params:
        outdir="{libname}/fastqc/",
        run_time = "02:09:00",
        cores="1",
        error_out_file = "error_files/fastqc.{libname}.txt",
        out_file = "stdout/fastqc.{libname}.txt",
    benchmark: "benchmarks/qc/fastqc.{libname}.{sample_label}.txt"
    shell:
        """
        module load fastqc;
        fastqc {input.fq1} --extract --outdir {params.outdir} -t {params.cores}
        fastqc {input.fq2} --extract --outdir {params.outdir} -t {params.cores}
        """

rule align_reads:
    input:
        fq1="{libname}/fastqs/ultraplex_demux_{sample_label}_Rev.Tr.fastq.gz",
        fq2="{libname}/fastqs/ultraplex_demux_{sample_label}_Fwd.Tr.fastq.gz",
    output:
        bam = "{libname}/bams/{sample_label}.Aligned.sortedByCoord.out.bam",
        bai = "{libname}/bams/{sample_label}.Aligned.sortedByCoord.out.bam.bai",
        unmapped1= "{libname}/bams/{sample_label}.Unmapped.out.mate1",
        unmapped2= "{libname}/bams/{sample_label}.Unmapped.out.mate2",
        log= "{libname}/bams/{sample_label}.Log.final.out",
    params:
        error_out_file = "error_files/align.{libname}.{sample_label}.err",
        out_file = "stdout/align.{libname}.{sample_label}.out",
        run_time = "02:00:00",
        memory = "40000",
        job_name = "align_reads",
        star_sjdb = config['STAR_DIR'],
        outprefix = "{libname}/bams/{sample_label}.",
        cores = "8",
    benchmark: "benchmarks/pre/align.{libname}.{sample_label}.txt"
    shell:        
        """
        module load star
        STAR \
            --alignEndsType Local \
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
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand zcat \
            --runMode alignReads \
            --runThreadN {params.cores}\

        module load samtools
        samtools index {output.bam}
        """

# rule align_reads_to_REPEAT:
#     input:
#         fq1="{libname}/fastqs/ultraplex_demux_{sample_label}_Rev.Tr.fastq.gz",
#         fq2="{libname}/fastqs/ultraplex_demux_{sample_label}_Fwd.Tr.fastq.gz",
#     output:
#         ubam = "{libname}/bams/repeat/{sample_label}.Aligned.out.bam",
#         unmapped1= "{libname}/bams/repeat/{sample_label}.Unmapped.out.mate1",
#         unmapped2= "{libname}/bams/repeat/{sample_label}.Unmapped.out.mate2",
#         log= "{libname}/bams/repeat/{sample_label}.Log.final.out",
#     params:
#         error_out_file = "error_files/{sample_label}_align_reads",
#         out_file = "stdout/{sample_label}_align_reads",
#         run_time = "06:40:00",
#         cores = "4",
#         memory = "10000",
#         job_name = "align_reads",
#         star_sjdb = config['STAR_REP'],
#         outprefix = "{libname}/bams/repeat/{sample_label}.",
#     benchmark: "benchmarks/align/{libname}.{sample_label}.align_reads.txt"
#     shell:
#         """
#         module load star ;
#         STAR \
#             --alignEndsType EndToEnd \
#             --genomeDir {params.star_sjdb} \
#             --genomeLoad NoSharedMemory \
#             --outBAMcompression 10 \
#             --outFileNamePrefix {params.outprefix} \
#             --outFilterMultimapNmax 30 \
#             --outFilterMultimapScoreRange 1 \
#             --outFilterScoreMin 10 \
#             --outFilterType BySJout \
#             --outReadsUnmapped Fastx \
#             --outSAMattrRGline ID:foo \
#             --outSAMattributes All \
#             --outSAMmode Full \
#             --outSAMtype BAM Unsorted \
#             --outSAMunmapped Within \
#             --outStd Log \
#             --readFilesIn {input.fq1} {input.fq2}\
#             --readFilesCommand zcat \
#             --runMode alignReads \
#             --runThreadN 8
#         """

# rule align_to_GENOME:
#     input:
#         fq1= "{libname}/bams/repeat/{sample_label}.Unmapped.out.mate1",
#         fq2= "{libname}/bams/repeat/{sample_label}.Unmapped.out.mate2",
#     output:
#         bam = "{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam",
#         bai = "{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam.bai",
#         unmapped1= "{libname}/bams/genome/{sample_label}.genome-mapped.Unmapped.out.mate1",
#         unmapped2= "{libname}/bams/genome/{sample_label}.genome-mapped.Unmapped.out.mate2",
#         log= "{libname}/bams/genome/{sample_label}.genome-mapped.Log.final.out",
#     params:
#         error_out_file = "error_files/{libname}.{sample_label}_align_reads_genome",
#         out_file = "stdout/{sample_label}_align_reads_genome",
#         run_time = "06:40:00",
#         cores = "4",
#         memory = "10000",
#         job_name = "align_reads",
#         star_sjdb = config['STAR_DIR'],
#         outprefix = "{libname}/bams/genome/{sample_label}.genome-mapped.",
#     benchmark: "benchmarks/align/{libname}.{sample_label}.align_reads.txt"
#     shell:
#         """
#         module load star ;
#         STAR \
#         --alignEndsType EndToEnd \
#         --genomeDir {params.star_sjdb} \
#         --genomeLoad NoSharedMemory \
#         --outBAMcompression 10 \
#         --outFileNamePrefix {params.outprefix} \
#         --outFilterMultimapNmax 1 \
#         --outFilterMultimapScoreRange 1 \
#         --outFilterScoreMin 10 \
#         --outFilterType BySJout \
#         --outReadsUnmapped Fastx \
#         --outSAMattrRGline ID:foo \
#         --outSAMattributes All \
#         --outSAMmode Full \
#         --outSAMtype BAM SortedByCoordinate \
#         --outSAMunmapped Within \
#         --outStd Log \
#         --readFilesIn {input.fq1} {input.fq2} \
#         --runMode alignReads \
#         --runThreadN 8

#         module load samtools
#         samtools index {output.bam}
#         """


rule umi_dedup:
    input:
        #bam="{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam",
        bam="{libname}/bams/{sample_label}.Aligned.sortedByCoord.out.bam",
        #bai="{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam.bai",
        bai="{libname}/bams/{sample_label}.Aligned.sortedByCoord.out.bam.bai"
    output:
        bam_dedup="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
    params:
        error_out_file = "error_files/dedup.{libname}.{sample_label}",
        out_file = "stdout/dedup.{libname}.{sample_label}",
        run_time = "06:40:00",
        cores = "4",
        memory = "10000",
        job_name = "sortbam",
        prefix='{libname}/bams/genome/{sample_label}.genome-mapped',
        java = config['JAVA_PATH'],
        umicollapse = config['UMICOLLAPSE_PATH']
    conda: 
        "envs/umi_tools.yaml"
    benchmark: "benchmarks/pre/umi_dedup.{libname}.{sample_label}.txt"
    shell:
        """
        {params.java} -server -Xms8G -Xmx8G -Xss20M -jar {params.umicollapse} bam -i {input.bam} -o {output.bam_dedup} --umi-sep : --two-pass --paired
        """
rule index_genome_bams:
    input:
        bam = "{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam"
    output:
        bai = "{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam.bai"
    params:
        error_out_file = "error_files/index_bam.{libname}.{sample_label}",
        out_file = "stdout/index_bam.{libname}.{sample_label}",
        run_time = "01:40:00",
        cores = "4",
        memory = "1000",
        job_name = "index_bam"
    benchmark: "benchmarks/pre/index_bam.{libname}.{sample_label}.txt"
    shell:
        "module load samtools;"
        "samtools index {input.bam};"
