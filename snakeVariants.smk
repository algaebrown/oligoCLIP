#snakemake -s snakeVariants.smk -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile config/preprocess_config/oligose_k562.yaml --use-conda --conda-prefix /tscc/nfs/home/hsher/snakeconda -np
import pandas as pd
workdir: config['WORKDIR']
SCRIPT_PATH = config['SCRIPT_PATH']
manifest = pd.read_table(config['MANIFEST'], index_col = False, sep = ',')
print(manifest)
barcode_df = pd.read_csv(config['barcode_csv'], header = None, sep = ':', names = ['barcode', 'RBP'])
# basic checking
assert not barcode_df['barcode'].duplicated().any()
assert not barcode_df['RBP'].duplicated().any() # cannot have any duplicated RBP names
assert not barcode_df['RBP'].str.contains(' ').any() # DO NOT CONTAIN white space lah
assert not manifest['fastq'].duplicated().any()
assert not manifest['libname'].str.contains(' ').any()
libnames = manifest['libname'].tolist() 

config['libnames'] = libnames
experiments = manifest['experiment'].tolist()
config['experiments'] = experiments
rbps = barcode_df['RBP'].tolist()
config['rbps'] = rbps
seq_df = '/tscc/nfs/home/hsher/scratch/ABC_DL/output/tsv/K562_rep6.RBFOX2.tsv' # need to be of the same as the feature annotatiions
rule all:
    input:
        expand("variants/{signal_type}/{libname}.{sample_label}.csv",
        signal_type = ['CITS'], 
        libname = ['K562_rep6'],
        sample_label = rbps,
        )
rule fetch_SNP:
    input:
        finemapped_windows = "DMM/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz"
    output:
        temp("variants/{signal_type}/{libname}.{sample_label}.{chr}.vcf")
    params:
        error_out_file = "error_files/fetch_snp.{signal_type}.{libname}.{sample_label}.{chr}",
        out_file = "stdout/fetch_snp.{signal_type}.{libname}.{sample_label}.{chr}",
        run_time = "3:20:00",
        cores = 1,
    container:
        "docker://miguelpmachado/bcftools:1.9-01"
    shell:
        """
        bcftools query -R {input.finemapped_windows} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' \
            https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.{wildcards.chr}.vcf.bgz > {output}
        """

rule fetch_sequence:
    input:
        subset_vcf="variants/{signal_type}/{libname}.{sample_label}.{chr}.vcf",
        seq_df=seq_df,
        feature_annotation=config['FEATURE_ANNOTATIONS'],
    output:
        temp("variants/{signal_type}/{libname}.{sample_label}.{chr}.csv")
    params:
        error_out_file = "error_files/fetch_sequence.{signal_type}.{libname}.{sample_label}.{chr}",
        out_file = "stdout/fetch_sequence.{signal_type}.{libname}.{sample_label}.{chr}",
        run_time = "03:20:00",
        cores = 1
    conda:
        "/tscc/nfs/home/hsher/projects/oligoCLIP/rules/envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/generate_variant_sequence.py {input.subset_vcf} {input.seq_df} {input.feature_annotation} \
            {output}
        """
rule combine_csv:
    input:
        expand("variants/{signal_type}/{libname}.{sample_label}.{chr}.csv",
            signal_type = ['{signal_type}'], 
            libname = ['{libname}'],
            sample_label = ['{sample_label}'],
            chr = [f'chr{i}' for i in list(range(1,23))+['X','Y']])
    output:
        "variants/{signal_type}/{libname}.{sample_label}.csv"
    params:
        error_out_file = "error_files/combine.{signal_type}.{libname}.{sample_label}",
        out_file = "stdout/combine.{signal_type}.{libname}.{sample_label}",
        run_time = "03:20:00",
        cores = 1
    shell:
        """
        awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} > {output}
        """