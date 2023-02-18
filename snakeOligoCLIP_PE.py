
from importlib.resources import path
import pandas as pd
import os
#snakemake -s snakeOligoCLIP_PE.py -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile config/preprocess_config/oligope_iter4.yaml --use-conda --conda-prefix /home/hsher/snakeconda -n
#snakemake -s snakeOligoCLIP_PE.py -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile config/preprocess_config/oligope_iter3_multimap.yaml --use-conda --conda-prefix /home/hsher/snakeconda -n  
#snakemake -s snakeOligoCLIP_PE.py -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" --configfile config/preprocess_config/oligope_sara_MN.yaml --use-conda --conda-prefix /home/hsher/snakeconda -n  
MANIFEST=config['MANIFEST']
SCRIPT_PATH=config['SCRIPT_PATH']
WORKDIR=config['WORKDIR']
UNINFORMATIVE_READ = 3 - int(config['INFORMATIVE_READ']) # whether read 1 or read 2 is informative
CHROM_SIZES = config['CHROM_SIZES']
workdir: WORKDIR
config['GENOME_FA'] = config['GENOMEFA']
R_EXE = config['R_EXE']

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

module DMN:
    snakefile:
        "rules/normalization_DMN.smk"
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

module analysis:
    snakefile:
        "rules/analysis.smk"
    config:
        config

module multimap:
    snakefile:
        "rules/multimap.smk"
    config:
        config

module merge_bw:
    snakefile:
        "rules/merge_bw.smk"
    config:
        config

module finemap:
    snakefile:
        "rules/finemap.smk"
    config:
        config

module clipper:
    snakefile:
        "rules/clipper.smk"
    config:
        config

module clipper_analysis:
    snakefile:
        "rules/clipper_analysis.smk"
    config:
        config



rule all:
    input:
        expand("{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam.bai", libname = libnames, sample_label = rbps),
        'QC/fastQC_basic_summary.csv',
        'QC/fastQC_passfail.csv',
        'QC/cutadapt_stat.csv',
        "QC/mapping_stats.csv",
        "QC/dup_level.csv",
        'QC/demux_read_count.txt',
        # expand("QC/nobarcode_blast_output/{libname}.blast.tsv", libname = libnames),
        # expand("QC/unmapped_blast_output/{libname}.{sample_label}.1.blast.tsv", libname = libnames, sample_label = rbps),
        # expand("QC/unmapped_blast_output/{libname}.{sample_label}.short.blast.tsv", libname = libnames, sample_label = rbps),
        expand("{libname}/bw/COV/{sample_label}.r1.{strand}.bw", libname = libnames, sample_label = rbps, strand = ['pos', 'neg']),
        expand("{libname}/bw_bg/COV/{sample_label}.r1.{strand}.bw", libname = libnames, sample_label = rbps, strand = ['pos', 'neg']),
        ############### skipper #####################
        # expand("output/enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_windows.tsv.gz",
        # libname = libnames,
        # clip_sample_label = list(set(rbps)-set([config['AS_INPUT']])), # cannot call on itself
        # bg_sample_label = [config['AS_INPUT']] if config['AS_INPUT'] else []
        # ,),
        # expand("internal_output/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
        # libname = libnames,
        # clip_sample_label = list(set(rbps)-set([config['AS_INPUT']])), # cannot call on itself
        # ),
        # expand("internal_output/enriched_re/{libname}.{sample_label}.enriched_re.tsv.gz",
        # libname = libnames,
        # sample_label = list(set(rbps)-set([config['AS_INPUT']])), # cannot call on itself
        # ),
        # expand("debug/counts/genome/tables/{libname}.{sample_label}.tsv.gz",
        # libname = libnames,
        # sample_label =['SLBP', 'PUM2', 'DDX3']
        # ),
        # expand("internal_output/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
        # libname = libnames,
        # sample_label = list(set(rbps)-set([config['AS_INPUT']])), 
        # signal_type = ['CITS', 'COV'] # cannot call on itself
        # ),
        # expand("internal_output/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
        # libname = libnames,
        # sample_label = config['RBP_TO_RUN_MOTIF'],
        # signal_type = ['CITS', 'COV']
        # ),
        ####### CLIPper ########
        # expand("output/CLIPper.{bg}/{libname}.{sample_label}.peaks.normed.compressed.annotate.bed",
        # bg = [config['AS_INPUT']] if config['AS_INPUT'] else [],
        # sample_label = list(set(rbps)-set([config['AS_INPUT']])),
        # libname = libnames
        # ),
        # expand("output/CLIPper-CC/{libname}.{sample_label}.peaks.normed.compressed.annotate.bed",
        # sample_label = list(set(rbps)-set([config['AS_INPUT']])),
        # libname = libnames
        # ),
        # expand("output/CLIPper-CC/{libname}.{sample_label}.peaks.normed.compressed.motif.svg",
        # sample_label = config['RBP_TO_RUN_MOTIF'],
        # libname = libnames
        # ),
        expand("internal_output/DMN/most_enriched_selected/{libname}.{sample_label}.enriched_window.tsv",
        libname = libnames,
        sample_label = list(set(rbps)-set([config['AS_INPUT']])),
        signal_type = ['CITS', 'COV'],
        strand = ['pos', 'neg']
        ),
        ####### DMN ########
        ### CC ###
        expand("internal_output/DMN/finemapping/mapped_sites_most_enriched/{signal_type}/{libname}.{sample_label}.finemapped_windows.{strand}.bw",
        libname = libnames,
        sample_label = list(set(rbps)-set([config['AS_INPUT']])),
        signal_type = ['CITS', 'COV'],
        strand = ['pos', 'neg']
        ),
        expand("internal_output/DMN/homer_most_enriched/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
        libname = libnames,
        sample_label = config['RBP_TO_RUN_MOTIF'],
        signal_type = ['CITS', 'COV']
        ),
        ### IgG ###
        expand("output/DMN/{bg_sample_label}/most_enriched_selected/{libname}.{clip_sample_label}.enriched_window.tsv",
        libname = libnames,
        clip_sample_label = list(set(rbps)-set([config['AS_INPUT']])),
        bg_sample_label = [config['AS_INPUT']] if config['AS_INPUT'] else []
        ),
        
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
use rule * from analysis as analysis_*

############# Quality control #################

use rule gather_trimming_stat from QC as qc_trim with:
    input:
        tr1=expand("QC/{libname}.Tr.metrics", libname = libnames)
    output:
        tr1='QC/cutadapt_stat.csv'

use rule gather_fastqc_report from QC as fastqc_gather with:
    input:
        expand("{libname}/fastqc/ultraplex_demux_{sample_label}_Rev.Tr_fastqc/fastqc_data.txt", libname = libnames, sample_label = rbps)+
        expand("{libname}/fastqc/ultraplex_demux_{sample_label}_Fwd.Tr_fastqc/fastqc_data.txt", libname = libnames, sample_label = rbps)

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

########## EXTRA DEBUGGING QC ############
use rule what_is_read_wo_barcode from QC
use rule blast_unmapped_reads from QC
use rule count_demultiplex_ultraplex from QC
use rule  blast_unmapped_reads_too_short from QC

### region caller ###
use rule * from normalization as skipper_*

### repeat caller ###
use rule * from repeat as re_*

############## DMN #################
use rule * from DMN as dmn_*

############## DMN #################
use rule * from clipper as clipper_*
use rule * from clipper_analysis

########## BIGWIGS ############
use rule extract_read_two from make_track as extract_r1 with: # oligoPE truncation is in read1
    input:
        bam="{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam"
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
########## MERGE BW ############
use rule * from merge_bw

########## DEBUG: FILTER MULTIMAPPED ############
use rule * from multimap

use rule get_nt_coverage from finemap with:
    input:
        windows = "internal_output/enriched_windows/{libname}.{sample_label}.enriched_windows.tsv.gz",
        clip_bw_pos = "{libname}/bw/{signal_type}/{sample_label}.r1.pos.bw",
        clip_bw_neg = "{libname}/bw/{signal_type}/{sample_label}.r1.neg.bw",
        input_bw_pos = "{libname}/bw_bg/{signal_type}/{sample_label}.r1.pos.bw",
        input_bw_neg = "{libname}/bw_bg/{signal_type}/{sample_label}.r1.neg.bw",

use rule get_nt_coverage_dmn from finemap with:
    input:
        windows = "internal_output/DMN/{libname}.{sample_label}.enriched_windows.tsv.gz",
        clip_bw_pos = "{libname}/bw/{signal_type}/{sample_label}.r1.pos.bw",
        clip_bw_neg = "{libname}/bw/{signal_type}/{sample_label}.r1.neg.bw",
        input_bw_pos = "{libname}/bw_bg/{signal_type}/{sample_label}.r1.pos.bw",
        input_bw_neg = "{libname}/bw_bg/{signal_type}/{sample_label}.r1.neg.bw",

use rule get_nt_coverage_dmn_most_enriched from finemap with:
    input:
        windows = "internal_output/DMN/most_enriched_selected/{libname}.{sample_label}.enriched_window.tsv",
        clip_bw_pos = "{libname}/bw/{signal_type}/{sample_label}.r1.pos.bw",
        clip_bw_neg = "{libname}/bw/{signal_type}/{sample_label}.r1.neg.bw",
        input_bw_pos = "{libname}/bw_bg/{signal_type}/{sample_label}.r1.pos.bw",
        input_bw_neg = "{libname}/bw_bg/{signal_type}/{sample_label}.r1.neg.bw",



use rule finemap_windows from finemap
use rule finemap_windows_dmn from finemap
use rule finemap_windows_dmn_most_enriched from finemap
use rule finemap_to_bedgraph from finemap