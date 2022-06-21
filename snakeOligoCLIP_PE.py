
from importlib.resources import path
import pandas as pd
import os

#snakemake -s snakeOligoCLIP_PE.py -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o /dev/null" --configfile config/preprocess_config/oligope.yaml --use-conda -n 

try:
    MANIFEST=config['MANIFEST']
    #SCRIPT_PATH=config['SCRIPT_PATH']
    WORKDIR=config['WORKDIR']
    workdir: WORKDIR

    print('WORKDIR', os.getcwd())

    manifest = pd.read_table(MANIFEST, index_col = False, sep = ',')
    barcode_df = pd.read_csv(config['barcode_csv'], header = None, sep = ':', names = ['barcode', 'RBP'])
    #print(manifest.head())

    # basic checking
    assert not barcode_df['barcode'].duplicated().any()
    assert not barcode_df['RBP'].duplicated().any() # cannot have any duplicated RBP names
    assert not barcode_df['RBP'].str.contains(' ').any() # DO NOT CONTAIN white space lah
    assert not manifest['fastq1'].duplicated().any()
    assert not manifest['fastq2'].duplicated().any()
    assert not manifest['libname'].str.contains(' ').any()

    libnames = manifest['libname'].tolist() 
    
    rbps = barcode_df['RBP'].tolist()
    # sample_labels = manifest.uid.tolist()
    # print(sample_labels)
    #sample_labels = ['Dan_singleplex_HEK293_rep1_RBFOX2','Dan_singleplex_HEK293_rep2_RBFOX2','ENCODE_Dan_singleplex_HEK293_rep1_RBFOX2','ENCODE_Dan_singleplex_HEK293_rep2_RBFOX2']
    
    try:
        os.mkdir('error_files')
    except:
        pass

    # all_rbfox = [s for s in sample_labels if 'RBFOX2' in s or '676' in s]
    # print(','.join(all_rbfox))
except Exception as e:
    print('config:', e)

module preprocess:
    snakefile:
        "rules/pe_preprocess.py"
    config: config

module QC:
    snakefile:
        "rules/QC.py"
    config:
        config



rule all:
    input:
        expand("{libname}/bams/genome/{sample_label}.genome-mapped.rmDup.Aligned.sortedByCoord.out.bam.bai", libname = libnames, sample_label = rbps),
        'QC/fastQC_basic_summary.csv',
        'QC/fastQC_passfail.csv',
        'QC/cutadapt_stat.csv',
        "QC/repeat_mapping_stats.csv",
        "QC/genome_mapping_stats.csv",
        "QC/dup_level.csv",
        expand("QC/nobarcode_blast_output/{libname}.blast.tsv", libname = libnames)

        
    output:
        "snakeCLIP.txt"
    params:
        error_out_file = "error_files/all",
        run_time = 1,
        cores = "1",
        memory = "20",
        job_name = "all"
    shell:
        """
        echo $(date)  > {output};
        echo created by HLH and the Yeo lab >> {output}
        """

use rule * from preprocess as pre_*

use rule gather_trimming_stat from QC as qc_trim with:
    input:
        tr1=expand("QC/{libname}.Tr.metrics", libname = libnames)
    output:
        tr1='QC/cutadapt_stat.csv'

use rule gather_fastqc_report from QC as fastqc_gather with:
    input:
        expand("{libname}/fastqc/ultraplex_demux_{sample_label}_Rev_fastqc/fastqc_data.txt", libname = libnames, sample_label = rbps)+
        expand("{libname}/fastqc/ultraplex_demux_{sample_label}_Fwd_fastqc/fastqc_data.txt", libname = libnames, sample_label = rbps)

use rule gather_mapstat from QC as mapstat_gather_repeat with:
    input:
        expand("{libname}/bams/repeat/{sample_label}.Log.final.out", libname = libnames, sample_label = rbps),
    output:
        "QC/repeat_mapping_stats.csv"

use rule gather_mapstat from QC as mapstat_gather_genome with:
    input:
        expand("{libname}/bams/genome/{sample_label}.genome-mapped.Log.final.out", libname = libnames, sample_label = rbps),
    output:
        "QC/genome_mapping_stats.csv"

use rule duplication_rate from QC as qc_duplication_rate with:
    input:
        dup=expand("{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam", libname = libnames, sample_label = rbps),
        rmdup=expand("{libname}/bams/genome/{sample_label}.genome-mapped.rmDup.Aligned.sortedByCoord.out.bam",libname = libnames, sample_label = rbps)
    output:
        "QC/dup_level.csv"

use rule what_is_read_wo_barcode from QC