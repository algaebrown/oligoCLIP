rule extract_multimap_and_uniquemap:
    input:
        "{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam"
    output:
        multi="debug/bam/{libname}.{sample_label}.multimap.bam",
        uniq="debug/bam/{libname}.{sample_label}.uniqmap.bam",
    params:
        run_time="1:00:00",
        error_out_file = "error_files/{libname}.{sample_label}.extractmap",
        out_file = "stdout/{libname}.{sample_label}.extractmap",
        cores = 1,
        memory = 10000,
    conda:
        "envs/bamtools"
    shell:
        """
        bamtools filter -in {input} -out {output.multi} -mapQuality "<=3"
        bamtools filter -in {input} -out {output.uniq} -mapQuality ">3"
        """

module skipper:
    snakefile:
        "skipper.smk"
    config:
        config

use rule partition_bam_reads from skipper as partition_map_by_status with:
    input:
        chrom_size = config['CHROM_SIZES'],
        bam = "debug/bam/{libname}.{sample_label}.{mapstat}.bam",        
        region_partition = config['PARTITION'],
    output:
        counts= "debug/counts/genome/vectors/{libname}.{sample_label}.{mapstat}.counts",
    params:
        error_out_file = "error_files/{libname}.{sample_label}.{mapstat}.partition_bam_reads.err",
        out_file = "stdout/{libname}.{sample_label}.{mapstat}.partition_bam_reads.out",
        run_time = "20:00",
        cores = "1",
        memory = 10000,
        replicate_label = "{libname}.{sample_label}.{mapstat}"
    benchmark: "benchmarks/counts/unassigned_experiment.{libname}.{sample_label}.{mapstat}.partition_bam_reads.txt"

# concat all the experiments of the same set into table
use rule make_genome_count_table from skipper as make_genome_count_bymap with:
    input:
        partition=config['PARTITION'],
        replicate_counts = lambda wildcards: expand("debug/counts/genome/vectors/{libname}.{sample_label}.{mapstat}.counts", 
            mapstat = ['multimap', 'uniqmap'], # TODO: make dictionary
            libname = [wildcards.libname],
            sample_label = [wildcards.sample_label]),
    output:
        count_table = "debug/counts/genome/tables/{libname}.{sample_label}.tsv.gz",
    params:
        error_out_file = "error_files/debug.{libname}.{sample_label}.make_count_table.err",
        out_file = "stdout/debug.{libname}.{sample_label}.make_count_table.out",
        run_time = "00:05:00",
        cores = "1",
        memory = 200,
    benchmark: "benchmarks/counts/debug.{libname}.{sample_label}.all_replicates.make_genome_count_table.txt"