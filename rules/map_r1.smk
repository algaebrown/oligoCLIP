# snakemake cannot have the same module and rule name even in submodules
module QC_1:
    snakefile:
        "QC.py"
    config:
        config

module make_track_1:
    snakefile:
        "rules/make_track.smk"
    config:
        config

rule align_to_GENOME_r1_only:
    input:
        fq1= "{libname}/bams/repeat/{sample_label}.Unmapped.out.mate1",
    output:
        bam = "{libname}/bams/genome_r1/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam",
        unmapped1= "{libname}/bams/genome_r1/{sample_label}.genome-mapped.Unmapped.out.mate1",
        log= "{libname}/bams/genome_r1/{sample_label}.genome-mapped.Log.final.out",
    params:
        error_out_file = "error_files/{libname}.{sample_label}_align_reads_genome",
        out_file = "stdout/{sample_label}_align_reads_genome",
        run_time = "06:40:00",
        cores = "4",
        memory = 10000,
        star_sjdb = config['STAR_DIR'],
        outprefix = "{libname}/bams/genome_r1/{sample_label}.genome-mapped.",
    benchmark: "benchmarks/align/{libname}.{sample_label}.align_reads.r1.txt"
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
        --readFilesIn {input.fq1} \
        --runMode alignReads \
        --runThreadN 8
        """

rule index_bam:
    input:
        "{anything}.bam"
    output:
        "{anything}.bam.bai"
    params:
        error_out_file = "error_files/{anything}_index_bam",
        out_file = "stdout/{anything}_index_bam",
        run_time = "40:00",
        cores = "1",
        memory = 10000,
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools index {input}
        """

rule umi_dedup:
    # TODO: containerize
    input:
        bam="{libname}/bams/genome_r1/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam",
        bai="{libname}/bams/genome_r1/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam.bai",
    output:
        bam_dedup="{libname}/bams/genome_r1/{sample_label}.genome-mapped.rmDup.Aligned.sortedByCoord.out.bam",
    params:
        error_out_file = "error_files/dedup.{libname}.{sample_label}",
        out_file = "stdout/{libname}.{sample_label}.index_reads",
        run_time = "06:40:00",
        cores = "4",
        memory = 10000,
        prefix='{libname}/bams/genome/{sample_label}.genome-mapped'
    container:
        "docker://howardxu520/skipper:umicollapse_1.0.0"
    benchmark: "benchmarks/align/dedup.{libname}.{sample_label}.txt"
    shell:
        """
        java -server -Xms8G -Xmx8G -Xss20M -jar /UMICollapse/umicollapse.jar bam -i {input.bam} -o {output.bam_dedup} --umi-sep : --two-pass
        """

######### mapstat ###########
use rule gather_mapstat from QC_1 as mapstat_gather_genome_r1 with:
    input:
        expand("{libname}/bams/genome_r1/{sample_label}.genome-mapped.Log.final.out", 
        libname = config['libnames'], 
        sample_label = config['rbps']),
    output:
        "QC/genome_r1_mapping_stats.csv"

######### tracks ###########
use rule COV_bam_to_bedgraph from make_track as COV_bedgraph_r1 with:
    input:
        bam="{libname}/bams/genome_r1/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam"
    output:
        pos="{libname}/bw/COV_r1/{sample_label}.r1.pos.bedgraph",
        neg="{libname}/bw/COV_r1/{sample_label}.r1.neg.bedgraph"