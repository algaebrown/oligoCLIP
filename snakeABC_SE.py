import pandas as pd

#snakemake -s snakeABC_SE.py -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile config/preprocess_config/oligose_single_slbp_k562.yaml --use-conda --conda-prefix /home/hsher/snakeconda -np
MANIFEST=config['MANIFEST']
SCRIPT_PATH=config['SCRIPT_PATH']
WORKDIR=config['WORKDIR']
UNINFORMATIVE_READ = 3 - int(config['INFORMATIVE_READ']) # whether read 1 or read 2 is informative
CHROM_SIZES = config['CHROM_SIZES']
workdir: WORKDIR
config['GENOME_FA'] = config['GENOMEFA']
R_EXE = config['R_EXE']
GFF_file='/home/hsher/gencode_coords/gencode.v38.primary_assembly.annotation.gff3'
DB_FILE='/home/hsher/scratch/gencode.v38.k562.ominclip.db'
GENOME_dir='/home/hsher/scratch/GRCh38.primary/'
GENOMEFA=config['GENOMEFA']

manifest = pd.read_table(MANIFEST, index_col = False, sep = ',')
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

if config['RBP_TO_RUN_MOTIF'] is None:
    config['RBP_TO_RUN_MOTIF'] = []

if len(rbps)==1:
    singleplex = True
else:
    singleplex = False
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

module preprocess:
    snakefile:
        "rules/se_preprocess.smk"
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

module DMN:
    snakefile:
        "rules/normalization_DMN.smk"
    config:
        config

module finemap:
    snakefile:
        "rules/finemap.smk"
    config:
        config

module repeat:
    snakefile:
        "rules/repeat.smk"
    config:
        config

module make_track:
    snakefile:
        config['MAKE_TRACK']
    config:
        config

module merge_bw:
    snakefile:
        "rules/merge_bw.smk"
    config:
        config

module analysis:
    snakefile:
        "rules/analysis.smk"
    config:
        config

module clipper:
    snakefile:
        "rules/clipper.smk"
    config:
        config

module compare:
    snakefile:
        "rules/compare.smk"
    config:
        config


def preprocess_outputs():
    ''' return preprocessing outputs'''
    outputs = expand("{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam.bai", libname = libnames, sample_label = rbps
    )+expand("{libname}/bw/COV/{sample_label}.{strand}.bw", libname = libnames, sample_label = rbps, strand = ['pos', 'neg']
    )+expand("{libname}/bw_bg/COV/{sample_label}.{strand}.bw", libname = libnames, sample_label = rbps, strand = ['pos', 'neg']
    )+['QC/fastQC_basic_summary.csv',
        'QC/fastQC_passfail.csv',
        'QC/cutadapt_stat.csv',
        "QC/mapping_stats.csv",
        "QC/dup_level.csv",
        'QC/demux_read_count.txt'
        ]+expand("output/counts/genome/vectors/{libname}.{sample_label}.counts",
        libname = libnames, sample_label = rbps)
    return outputs
def skipper_outputs():
    # skipper
    outputs = expand("internal_output/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
    libname = libnames,
    clip_sample_label = list(set(rbps)-set([config['AS_INPUT']])), # cannot call on itself
    )+expand("internal_output/enriched_windows/{libname}.{clip_sample_label}.dlogs.enriched_windows.tsv.gz",
    libname = libnames,
    clip_sample_label = list(set(rbps)-set([config['AS_INPUT']])), # cannot call on itself
    )+expand("internal_output/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
    libname = libnames,
    sample_label = list(set(rbps)-set([config['AS_INPUT']])), 
    signal_type = ['CITS', 'COV'] # cannot call on itself
    )+expand("internal_output/enriched_re/{libname}.{sample_label}.enriched_re.tsv.gz",
    libname = libnames,
    sample_label = list(set(rbps)-set([config['AS_INPUT']])), # cannot call on itself
    )+expand("internal_output/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
    libname = libnames,
    sample_label = config['RBP_TO_RUN_MOTIF'],
    signal_type = ['CITS', 'COV']
    )
    return outputs

def DMN_outputs():
    outputs = expand("internal_output/DMN/{libname}.{clip_sample_label}.mixture_weight.tsv",
    libname = libnames,
    clip_sample_label = list(set(rbps)-set([config['AS_INPUT']])), # cannot call on itself
    )+expand("internal_output/DMN/most_enriched_selected/{libname}.{sample_label}.enriched_window.tsv",
    sample_label = list(set(rbps)-set([config['AS_INPUT']])),
    libname = libnames
    )+expand("internal_output/DMN/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.{strand}.bw",
    libname = libnames,
    sample_label = list(set(rbps)-set([config['AS_INPUT']])),
    signal_type = ['CITS', 'COV'],
    strand = ['pos', 'neg']
    )+expand("internal_output/DMN/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
    libname = libnames,
    sample_label = config['RBP_TO_RUN_MOTIF'],
    signal_type = ['CITS', 'COV']
    )
    return output
def clipper_outputs():
    outputs = expand("output/CLIPper.{bg}/{libname}.{sample_label}.peaks.normed.compressed.bed",
        bg = [config['AS_INPUT']] if config['AS_INPUT'] else [],
        sample_label = list(set(rbps)-set([config['AS_INPUT']])),
        libname = libnames
        )+expand("output/CLIPper-CC/{libname}.{sample_label}.peaks.normed.compressed.bed",
        sample_label = list(set(rbps)-set([config['AS_INPUT']])),
        libname = libnames
        )
    return outputs

def comparison_outputs():
    outputs = expand("comparison/piranha/CC/{libname}.{sample_label}.bed",
        libname = libnames,
        sample_label =list(set(rbps)-set([config['AS_INPUT']])),
    )+expand("comparison/omniCLIP/output/{libname}.{sample_label}.omniclip_done.txt",
        libname = libnames,
        sample_label = list(set(rbps)-set([config['AS_INPUT']]))
    )+expand("comparison/pureclip/{libname}.{sample_label}.bind.bed",
        libname = libnames,
        sample_label = list(set(rbps)-set([config['AS_INPUT']]))
    )
    return outputs
def get_output(singleplex, clipper, skipper, comparison):
    output = preprocess_outputs()
    if singleplex:
        return output
    output += DMN_outputs()
    if skipper:
        output += skipper_outputs()
    if clipper:
        output += clipper_outputs()
    if comparison:
        output += comparison_outputs()
    return output

rule all:
    input:
        get_output(singleplex, config['run_clipper'], config['run_skipper'], config['run_comparison'])
    params:
        error_out_file = "error_files/all",
        run_time = "00:04:00",
        cores = "1",
        memory = "20",
        job_name = "all",
        out_file = "all"
    shell:
        "echo $(date)  > {output};"

############## PREPROCESS #################
use rule * from preprocess as pre_*

############## QUALITY CONTROL #################

use rule gather_trimming_stat from QC as qc_trim with:
    input:
        tr1=expand("QC/{libname}.Tr.metrics", libname = libnames)
    output:
        tr1='QC/cutadapt_stat.csv'

use rule gather_fastqc_report from QC as fastqc_gather with:
    input:
        expand("{libname}/fastqc/ultraplex_demux_{sample_label}.rev_fastqc/fastqc_data.txt", 
        libname = libnames, 
        sample_label = rbps)

use rule gather_mapstat from QC as mapstat_gather_repeat with:
    input:
        expand("{libname}/bams/{sample_label}.Log.final.out", libname = libnames, sample_label = rbps),
    output:
        "QC/mapping_stats.csv"

use rule duplication_rate from QC as qc_duplication_rate with:
    input:
        dup=expand("{libname}/bams/{sample_label}.Aligned.sortedByCoord.out.bam", libname = libnames, sample_label = rbps),
        rmdup=expand("{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",libname = libnames, sample_label = rbps)
    output:
        "QC/dup_level.csv"

use rule count_demultiplex_ultraplex from QC with:
    input:
        fq1=expand("{libname}/fastqs/ultraplex_demux_{sample_label}.fastq.gz", 
        libname = libnames, sample_label = rbps)
    

############## SKIPPER: GENOME #################
use rule * from normalization as skipper_*
use rule * from finemap as fine_*
############## SKIPPER: REPEAT #################
use rule * from repeat as re_*

############## DMN #################
use rule * from DMN as dmn_*

############## DMN #################
use rule * from clipper as clipper_*

############## BIGWIGS #################
use rule CITS_bam_to_bedgrah from make_track as CITS_bedgraph with:
    input:
        bam="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam"
    output:
        pos="{libname}/bw/CITS/{sample_label}.pos.bedgraph",
        neg="{libname}/bw/CITS/{sample_label}.neg.bedgraph"
use rule COV_bam_to_bedgraph from make_track as COV_bedgraph with:
    input:
        bam="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam"
    output:
        pos="{libname}/bw/COV/{sample_label}.pos.bedgraph",
        neg="{libname}/bw/COV/{sample_label}.neg.bedgraph"

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

########## MERGE BW ############
use rule merge_other_bw_as_bg from merge_bw with:
    input:
        bws=lambda wildcards: expand("{libname}/bw/{signal_type}/{sample_label}.{strand}.bw",
            libname = [wildcards.libname],
            sample_label = list(set(rbps)-set([wildcards.sample_label])-set([config['AS_INPUT']])),
            signal_type = [wildcards.signal_type],
            strand = [wildcards.strand]
        )
    output:
        tem = "{libname}/bw_bg/{signal_type}/{sample_label}.{strand}.temp.bedgraph",
        merged_bedgraph = "{libname}/bw_bg/{signal_type}/{sample_label}.{strand}.bedgraph",

########## HOMER ############

use rule * from analysis


######## COMPARE ###########
use rule * from compare