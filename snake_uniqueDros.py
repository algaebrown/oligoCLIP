# find unique drospophila mapped reads
#snakemake -s snake_uniqueDros.py -j 16 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --directory=/home/hsher/scratch/katie_drosphila --keep-going -n

sample_labels,= glob_wildcards("bams/genome/{sample_label}.genome-mapped.Aligned.out.bam")

rule all:
    input:
        expand('readcount/{sample_label}.count', sample_label = sample_labels),
    output:
        'QC/dros_uniquestat.txt'
    params:
        error_out_file = "error_files/mapstat",
        run_time = "00:40:00",
        cores = "1",
        memory = "10000",
        job_name = "gather_stat",
    shell:
        """
        awk -v OFS=';' '{{print FILENAME, $0}}' {input} > {output}
        """
    
    
rule find_mapped_read_id:
    input:
        genome="bams/genome/{sample_label}.genome-mapped.Aligned.out.bam",
        repeat="bams/repeat/{sample_label}.Aligned.out.bam",
        dros="bams/dros/{sample_label}.Aligned.out.bam"
    output:
        dros='readid/{sample_label}.dros.id',
        genome='readid/{sample_label}.genome.id',
        repeat='readid/{sample_label}.repeat.id'
    params:
        error_out_file = "error_files/mapstat",
        run_time = "00:40:00",
        cores = "1",
        memory = "10000",
        job_name = "gather_stat",
    shell:
        """
        module load samtools
        samtools view -F 4 {input.genome} | cut -f 1 | sort > {output.genome}
        samtools view -F 4 {input.dros} | cut -f 1 | sort > {output.dros}
        samtools view -F 4 {input.repeat} | cut -f 1 | sort > {output.repeat}
        """

rule find_unique_number:
    input:
        dros='readid/{sample_label}.dros.id',
        genome='readid/{sample_label}.genome.id',
        repeat='readid/{sample_label}.repeat.id'
    output:
        'readcount/{sample_label}.count'
    params:
        error_out_file = "error_files/mapstat",
        run_time = "00:40:00",
        cores = "1",
        memory = "10000",
        job_name = "gather_stat",
    shell:
        """
        diff  {input.dros}  {input.genome} | grep '<' | wc -l > {output}
        """