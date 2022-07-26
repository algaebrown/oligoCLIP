# rules to calculate variability of n replicates
# input: 
# - bed files of n IPs
# - bam files of n IPs (to count read)

#print(config['data'])
workdir: config['WORKDIR']

def get_all_sample(config = config):
    return [s['sample_label'] for s in config['data']]
def get_beds(sample_label, config = config):
    ''' return beds of replicates, given sample_label '''
    the_sample = [s for s in config['data'] if s['sample_label']==sample_label]
    print(the_sample)
    assert len(the_sample)==1
    the_sample = the_sample[0]
    return [r['bed'] for r in the_sample['reps']]
def get_bams(sample_label, config = config):
    ''' return bam of replicates, given sample_label '''
    the_sample = [s for s in config['data'] if s['sample_label']==sample_label]
    print(the_sample)
    assert len(the_sample)==1
    the_sample = the_sample[0]
    return [r['bam'] for r in the_sample['reps']]

rule all:
    input:
        expand("merged_bed/{sample_label}.bed", sample_label = get_all_sample()),
        expand("read_count/{sample_label}.bed", sample_label = get_all_sample()),
    output:
        "var.txt"
    params:
        cores="1",
        run_time = "01:00:00",
    shell:
        """
        cat DONE > {output}
        """

rule merge_bed_between_reps:
    input:
        beds = lambda wildcards: get_beds(wildcards.sample_label)
    output:
        merged_beds = "merged_bed/{sample_label}.bed"
    params:
        cores="1",
        run_time = "01:00:00",
    shell:
        """
        module load bedtools
        cat {input.beds} | bedtools sort | bedtools merge -s -c 4,5,6 -o collapse,collapse,distinct > {output}
        """
rule count_reads:
    input:
        merged_beds = "merged_bed/{sample_label}.bed",
        bams = lambda wildcards: get_bams(wildcards.sample_label)
    output:
        counts = "read_count/{sample_label}.bed"
    params:
        cores="1",
        run_time = "01:00:00",
    shell:
        """
        module load bedtools
        bedtools multicov -bams {input.bams} -bed {input.merged_beds} -s > {output.counts}
        """
