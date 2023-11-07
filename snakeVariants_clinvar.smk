#snakemake -s snakeVariants_clinvar.smk -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile config/preprocess_config/oligose_k562.yaml --use-conda --conda-prefix /home/hsher/snakeconda -np
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
seq_df = '/home/hsher/scratch/ABC_DL/output/tsv/K562_rep6.RBFOX2.tsv' # need to be of the same as the feature annotatiions
VCF='/projects/ps-yeolab5/hsher/clinvar/clinvar.vcf.gz'
rule all:
    input:
        expand("variants_clinvar/{signal_type}/{libname}.{sample_label}.csv",
        signal_type = ['CITS'], 
        libname = ['K562_rep6'],
        sample_label = ['RBFOX2'],
        )

rule reannotate:
    input:
        vcf=VCF,
        rename="/home/hsher/projects/oligoCLIP/utils/rename_chr.txt"
    output:
        VCF.replace('.vcf.gz', '.rename.vcf.gz')
    params:
        error_out_file = "error_files/rename_chr",
        out_file = "stdout/rename_chr",
        run_time = "3:20:00",
        cores = 1,
    shell:
        """
        module load bcftools
        bcftools annotate --rename-chrs {input.rename} \
            {input.vcf} -Oz -o {output}
        bcftools index {output}
        """
rule fetch_SNP:
    input:
        finemapped_windows = "DMM/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
        vcf = VCF.replace('.vcf.gz', '.rename.vcf.gz')
    output:
        "variants_clinvar/{signal_type}/{libname}.{sample_label}.vcf"
    params:
        error_out_file = "error_files/fetch_snp.{signal_type}.{libname}.{sample_label}",
        out_file = "stdout/fetch_snp.{signal_type}.{libname}.{sample_label}",
        run_time = "3:20:00",
        cores = 1,
    shell:
        """
        module load bcftools 
        bcftools query -R {input.finemapped_windows} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNDN\t%INFO/CLNVC\t%INFO/CLNSIG\n' \
            {input.vcf} > {output}
        """

rule fetch_sequence:
    input:
        subset_vcf="variants_clinvar/{signal_type}/{libname}.{sample_label}.vcf",
        seq_df=seq_df,
        feature_annotation=config['FEATURE_ANNOTATIONS'],
    output:
        "variants_clinvar/{signal_type}/{libname}.{sample_label}.csv"
    params:
        error_out_file = "error_files/fetch_sequence.{signal_type}.{libname}.{sample_label}",
        out_file = "stdout/fetch_sequence.{signal_type}.{libname}.{sample_label}",
        run_time = "01:20:00",
        cores = 1
    conda:
        "/home/hsher/projects/oligoCLIP/rules/envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/generate_variant_sequence_clinvar.py {input.subset_vcf} {input.seq_df} {input.feature_annotation} \
            {output}
        """
