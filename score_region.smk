#snakemake -s score_region.smk -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile config/preprocess_config/oligose_k562.yaml --use-conda --conda-prefix /tscc/nfs/home/hsher/snakeconda -np

import pandas as pd
from pathlib import Path
workdir: '/tscc/nfs/home/hsher/scratch/count_A3SS'
rootdir = Path('/tscc/nfs/home/hsher/scratch/ABC_2rep')
#cat /tscc/nfs/home/hsher/scratch/encode_kd/output/rMATs/SF3B4_K562/regions/A3SS/*/long_ss3.bed > ~/scratch/long_ss3.bed
#cat /tscc/nfs/home/hsher/scratch/encode_kd/output/rMATs/SF3B4_K562/regions/A3SS/*/long_5intron.bed > ~/scratch/long_5intron.bed
# cat /tscc/nfs/home/hsher/scratch/encode_kd/output/rMATs/PRPF8_K562/regions/A5SS/*/long_3intron.bed > ~/scratch/A5SS_long_3intron.bed
# cat /tscc/nfs/home/hsher/scratch/encode_kd/output/rMATs/PRPF8_K562/regions/A5SS/*/exon_diff_interval.bed > ~/scratch/A5SS_exon_diff_interval.bed
# cat /tscc/nfs/home/hsher/scratch/encode_kd/output/rMATs/RBFOX2_K562/regions/SE/*/casette_exon.bed > ~/scratch/casette_exon.bed
# cat /tscc/nfs/home/hsher/scratch/encode_kd/output/rMATs/RBFOX2_K562/regions/SE/*/casette_5intron.bed > ~/scratch/casette_5intron.bed
# cat /tscc/nfs/home/hsher/scratch/encode_kd/output/rMATs/RBFOX2_K562/regions/SE/*/casette_3intron.bed > ~/scratch/casette_3intron.bed
# cat /tscc/nfs/home/hsher/scratch/encode_kd/output/rMATs/RBFOX2_K562/regions/SE/*/exon5.bed > ~/scratch/exon5.bed
# cat /tscc/nfs/home/hsher/scratch/encode_kd/output/rMATs/RBFOX2_K562/regions/SE/*/exon3.bed > ~/scratch/exon3.bed



REGIONS=['/tscc/nfs/home/hsher/scratch/long_ss3.bed',
'/tscc/nfs/home/hsher/scratch/long_5intron.bed',
'/tscc/nfs/home/hsher/scratch/exon_diff_interval.bed',
'/tscc/nfs/home/hsher/scratch/A5SS_long_3intron.bed',
'/tscc/nfs/home/hsher/scratch/A5SS_exon_diff_interval.bed',
'/tscc/nfs/home/hsher/scratch/exon3.bed',
'/tscc/nfs/home/hsher/scratch/exon5.bed',
'/tscc/nfs/home/hsher/scratch/casette_exon.bed',
'/tscc/nfs/home/hsher/scratch/casette_5intron.bed',
'/tscc/nfs/home/hsher/scratch/casette_3intron.bed'
]
REGION_NAMES = [Path(i).name.replace('.bed','') for i in REGIONS]

REGIONS_DICT = dict(zip(REGION_NAMES, REGIONS))

MANIFEST=config['MANIFEST']
SCRIPT_PATH=config['SCRIPT_PATH']
UNINFORMATIVE_READ = 3 - int(config['INFORMATIVE_READ']) # whether read 1 or read 2 is informative
CHROM_SIZES = config['CHROM_SIZES']
DB_FILE=config['DB_FILE']
GENOME_dir=config['GENOME_dir']
GENOMEFA=config['GENOMEFA']

manifest = pd.read_table(MANIFEST, index_col = False, sep = ',')
print(manifest)
barcode_df = pd.read_csv(config['barcode_csv'], header = None, sep = ':', names = ['barcode', 'RBP'])
# basic checking
assert not barcode_df['barcode'].duplicated().any()
assert not barcode_df['RBP'].duplicated().any() # cannot have any duplicated RBP names
assert not barcode_df['RBP'].str.contains(' ').any() # DO NOT CONTAIN white space lah
assert not manifest['fastq'].duplicated().any()
assert not manifest['libname'].str.contains(' ').any()
libnames = manifest['libname'].tolist() 

config['libnames'] = ['K562_rep6']#libnames
experiments = manifest['experiment'].tolist()
config['experiments'] = experiments
rbps = barcode_df['RBP'].tolist()
config['rbps'] = rbps

print(f'RBPs: {rbps}',
    f'experiments:{experiments}',
    f'libnames:{libnames}')

# module normalization:
#     snakefile:
#         "rules/skipper.smk"
#     config:
#         config

# module DMN:
#     snakefile:
#         "rules/normalization_DMN.smk"
#     config:
#         config

rule all:
    input:
        expand("{region_name}/counts/genome/megatables/{libname}.tsv.gz",
        region_name = REGION_NAMES,
        libname = libnames)

rule partition_bam_reads:
    input:
        chrom_size = config['CHROM_SIZES'],
        bam = rootdir/"{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",        
        region_partition = lambda wildcards: REGIONS_DICT[wildcards.region_name]
    output:
        counts= "{region_name}/counts/genome/vectors/{libname}.{sample_label}.counts",
    params:
        error_out_file = "error_files/partition_bam_reads.{libname}.{sample_label}.{region_name}.err",
        out_file = "stdout/partition_bam_reads.{libname}.{sample_label}.{region_name}.out",
        run_time = "20:00",
        cores = "1",
        memory = "10000",
        replicate_label = "{libname}.{sample_label}",
        uninformative_read = config['UNINFORMATIVE_READ']
    benchmark: "benchmarks/counts/unassigned_experiment.{libname}.{sample_label}.{region_name}.partition_bam_reads.txt"
    shell:
        "bedtools bamtobed -i {input.bam} | awk '($1 != \"chrEBV\") && ($4 !~ \"/{params.uninformative_read}$\")' | bedtools flank -s -l 1 -r 0 -g {input.chrom_size} -i - | bedtools shift -p -1 -m 1 -g {input.chrom_size} -i - | bedtools coverage -counts -s -a {input.region_partition} -b - | cut -f 7 | awk 'BEGIN {{print \"{params.replicate_label}\"}} {{print}}' > {output.counts};"

rule make_window_by_barcode_table:
    input:
        counts = expand("{region_name}/counts/genome/vectors/{libname}.{sample_label}.counts",
            libname = ["{libname}"],
            sample_label = rbps,
            region_name = ["{region_name}"]),
    output:
        counts = "{region_name}/counts/genome/megatables/{libname}.tsv.gz",
    params:
        error_out_file = "error_files/window_by_barcode_table.{libname}.{region_name}.err",
        out_file = "stdout/window_by_barcode_table.{libname}.{region_name}.out",
        run_time = "20:00",
        cores = 1
    shell:
        """
        paste -d '\t' {input.counts} | gzip  > {output.counts}
        """