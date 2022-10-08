
from importlib.resources import path
import pandas as pd
import os

#snakemake -s snakeOligoCLIP_PE.py -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile config/preprocess_config/oligope_iter2_tile.yaml --use-conda  --conda-prefix /home/hsher/scratch/oligo_PE/.snakemake/conda -n  
#snakemake --cluster "sbatch -t {params.run_time} -p home-yeo -e {params.error_out_file} -o {params.out_file} -T {params.cores}"
#snakemake -s snakeOligoCLIP_PE.py -j 12 --cluster "sbatch -t {params.run_time} -p home-yeo -e {params.error_out_file} -o {params.out_file} -T {params.cores}" --configfile config/preprocess_config/oligope_iter2_simulate.yaml --use-conda  --conda-prefix /home/hsher/scratch/oligo_PE/.snakemake/conda -n  
#snakemake -s snakeOligoCLIP_PE.py -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile config/preprocess_config/oligope_sara_CN.yaml --use-conda --conda-prefix /home/hsher/snakeconda -n  
MANIFEST=config['MANIFEST']
SCRIPT_PATH=config['SCRIPT_PATH']
WORKDIR=config['WORKDIR']
workdir: WORKDIR

manifest = pd.read_table(MANIFEST, index_col = False, sep = ',')
barcode_df = pd.read_csv(config['barcode_csv'], header = None, sep = ':', names = ['barcode', 'RBP'])

# basic checking
assert not barcode_df['barcode'].duplicated().any()
assert not barcode_df['RBP'].duplicated().any() # cannot have any duplicated RBP names
assert not barcode_df['RBP'].str.contains(' ').any() # DO NOT CONTAIN white space lah
assert not manifest['fastq1'].duplicated().any()
assert not manifest['fastq2'].duplicated().any()
assert not manifest['libname'].str.contains(' ').any()

libnames = manifest['libname'].tolist() 

config['libnames'] = libnames
experiments = manifest['experiment'].tolist()
config['experiments'] = experiments
rbps = barcode_df['RBP'].tolist()
config['rbps'] = rbps
# sample_labels = manifest.uid.tolist()
# print(sample_labels)
#sample_labels = ['Dan_singleplex_HEK293_rep1_RBFOX2','Dan_singleplex_HEK293_rep2_RBFOX2','ENCODE_Dan_singleplex_HEK293_rep1_RBFOX2','ENCODE_Dan_singleplex_HEK293_rep2_RBFOX2']

# making the error files directory
try:
    os.mkdir('error_files')
except:
    pass

# making the stdout directory
try:
    os.mkdir('stdout')
except:
    pass


config['GENOME_FA'] = config['GENOMEFA']
R_EXE = config['R_EXE']
# all_rbfox = [s for s in sample_labels if 'RBFOX2' in s or '676' in s]
# print(','.join(all_rbfox))

module preprocess:
    snakefile:
        "rules/pe_preprocess.py"
    config: config

module mapr1:
    snakefile:
        "rules/map_r1.py"
    config: config

module QC:
    snakefile:
        "rules/QC.py"
    config:
        config

module normalization:
    snakefile:
        "rules/normalization.py"
    config:
        config

module make_track:
    snakefile:
        config['MAKE_TRACK']
    config:
        config

module analysis:
    snakefile:
        "rules/analysis.py"
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
        'QC/demux_read_count.txt',
        expand("QC/nobarcode_blast_output/{libname}.blast.tsv", libname = libnames),
        expand("QC/unmapped_blast_output/{libname}.{sample_label}.1.blast.tsv", libname = libnames, sample_label = rbps),
        expand("QC/unmapped_blast_output/{libname}.{sample_label}.short.blast.tsv", libname = libnames, sample_label = rbps),
        expand("output/enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_windows.tsv.gz",
        libname = libnames,
        clip_sample_label = list(set(rbps)-set([config['AS_INPUT']])), # cannot call on itself
        bg_sample_label = [config['AS_INPUT']] if config['AS_INPUT'] else []
        ,),
        expand("internal_output/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
        libname = libnames,
        clip_sample_label = list(set(rbps)-set([config['AS_INPUT']])), # cannot call on itself
        ),
        expand("{libname}/bw/{type}/{sample_label}.r1.{strand}.bw", libname = libnames, sample_label = rbps, strand = ['pos', 'neg'], type = ['CITS', 'COV', 'COV_r1']),
        "QC/genome_r1_mapping_stats.csv",
        expand("output/homer/results/{libname}.{clip_sample_label}.{region}/homerMotifs.all.motifs",
        libname = libnames,
        clip_sample_label = list(set(rbps)-set([config['AS_INPUT']])),
        region = ['CDS', 'UTR5', 'INTRON', 'UTR3']
        )

        
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
use rule * from mapr1 as r1_*
use rule * from analysis as analysis_*

############# Quality control #################

use rule gather_trimming_stat from QC as qc_trim with:
    input:
        tr1=expand("QC/{libname}.Tr.metrics", libname = libnames)
    output:
        tr1='QC/cutadapt_stat.csv'

# use rule gather_trimming_stat from QC as qc_barcode_trim with:
#     input:
#         tr1=expand("QC/{libname}.{sample_label}.barcodeTr.metrics", 
#         libname = libnames,
#         sample_label = rbps)
#     output:
#         tr1=expand('QC/cutadapt_barcode.csv')

use rule gather_fastqc_report from QC as fastqc_gather with:
    input:
        expand("{libname}/fastqc/ultraplex_demux_{sample_label}_Rev.Tr_fastqc/fastqc_data.txt", libname = libnames, sample_label = rbps)+
        expand("{libname}/fastqc/ultraplex_demux_{sample_label}_Fwd.Tr_fastqc/fastqc_data.txt", libname = libnames, sample_label = rbps)

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
use rule blast_unmapped_reads from QC
use rule count_demultiplex_ultraplex from QC
use rule  blast_unmapped_reads_too_short from QC

### region caller ###
use rule * from normalization as skipper_*

### make CITS and COVERAGE track
use rule extract_read_two from make_track as extract_r1 with: # oligoPE truncation is in read1
    input:
        bam="{libname}/bams/genome/{sample_label}.genome-mapped.rmDup.Aligned.sortedByCoord.out.bam"
    output:
        read2="{libname}/bw/{sample_label}.r2.bam",
        read1="{libname}/bw/{sample_label}.r1.bam"

use rule strand_specific_bam from make_track as strand_specific_bam with:
    input:
        bam = "{libname}/bw/{sample_label}.r1.bam"
    output:
        pos_bam = "{libname}/bw/{sample_label}.r1.pos.bam",
        neg_bam = "{libname}/bw/{sample_label}.r1.neg.bam"

use rule CITS_bam_to_bedgrah from make_track as CITS_bedgraph with:
    input:
        bam="{libname}/bw/{sample_label}.r1.bam"
    output:
        pos="{libname}/bw/CITS/{sample_label}.r1.pos.bedgraph",
        neg="{libname}/bw/CITS/{sample_label}.r1.neg.bedgraph"
use rule COV_bam_to_bedgraph from make_track as COV_bedgraph with:
    input:
        bam="{libname}/bw/{sample_label}.r1.bam"
    output:
        pos="{libname}/bw/COV/{sample_label}.r1.pos.bedgraph",
        neg="{libname}/bw/COV/{sample_label}.r1.neg.bedgraph"

rule bedgraph_to_bw:
    input:
        bedgraph="{something}.bedgraph",
    output:
        bw="{something}.bw"
    params:
        run_time="6:00:00",
        chr_size=config['CHROM_SIZES'],
        error_out_file = "error_files/{something}.bedgraph_to_bw",
        out_file = "stdout/{something}.bedgraph_to_bw",
        cores = 1,
    shell:
        """
        module load ucsctools
        bedGraphToBigWig {input.bedgraph} {params.chr_size} {output.bw}
        """

