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
        bai = "{libname}/bams/genome_r1/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam.bai",
        unmapped1= "{libname}/bams/genome_r1/{sample_label}.genome-mapped.Unmapped.out.mate1",
        log= "{libname}/bams/genome_r1/{sample_label}.genome-mapped.Log.final.out",
    params:
        error_out_file = "error_files/{libname}.{sample_label}_align_reads_genome",
        out_file = "stdout/{sample_label}_align_reads_genome",
        run_time = "06:40:00",
        cores = "4",
        memory = "10000",
        job_name = "align_reads",
        star_sjdb = config['STAR_DIR'],
        outprefix = "{libname}/bams/genome_r1/{sample_label}.genome-mapped.",
    benchmark: "benchmarks/align/{libname}.{sample_label}.align_reads.r1.txt"
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
        --readFilesIn {input.fq1} \
        --runMode alignReads \
        --runThreadN 8

        module load samtools
        samtools index {output.bam}
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
        {params.java} -server -Xms8G -Xmx8G -Xss20M -jar {params.umicollapse} bam -i {input.bam} -o {output.bam_dedup} --umi-sep : --two-pass
        """
rule index_genome_bams:
    input:
        bam = "{libname}/bams/genome_r1/{sample_label}.genome-mapped.rmDup.Aligned.sortedByCoord.out.bam"
    output:
        bai = "{libname}/bams/genome_r1/{sample_label}.genome-mapped.rmDup.Aligned.sortedByCoord.out.bam.bai"
    params:
        error_out_file = "error_files/index_bam.{libname}.{sample_label}",
        out_file = "stdout/index_bam.{libname}.{sample_label}",
        run_time = "01:40:00",
        cores = "4",
        memory = "1000",
        job_name = "index_bam"
    benchmark: "benchmarks/align/{libname}.{sample_label}.index_bam.txt"
    shell:
        "module load samtools;"
        "samtools index {input.bam};"

######### mapstat ###########
use rule gather_mapstat from QC_1 as mapstat_gather_genome_r1 with:
    input:
        expand("{libname}/bams/genome_r1/{sample_label}.genome-mapped.Log.final.out", 
        libname = config['libnames'], 
        sample_label = config['rbps']),
    output:
        "QC/genome_r1_mapping_stats.csv"

######### tracks ###########
use rule COV_bam_to_bedgraph from make_track_1 as COV_bedgraph_r1 with:
    input:
        bam="{libname}/bams/genome_r1/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam"
    output:
        pos="{libname}/bw/COV_r1/{sample_label}.r1.pos.bedgraph",
        neg="{libname}/bw/COV_r1/{sample_label}.r1.neg.bedgraph"