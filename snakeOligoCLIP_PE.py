
from importlib.resources import path
import pandas as pd
import os

#snakemake -s snakeOligoCLIP_PE.py -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o /dev/null" --configfile config/preprocess_config/oligope.yaml --use-conda  --conda-prefix /home/hsher/scratch/oligo_PE/.snakemake/conda -n 

try:
    MANIFEST=config['MANIFEST']
    SCRIPT_PATH=config['SCRIPT_PATH']
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
    experiments = manifest['experiment'].tolist()
    
    rbps = barcode_df['RBP'].tolist()
    # sample_labels = manifest.uid.tolist()
    # print(sample_labels)
    #sample_labels = ['Dan_singleplex_HEK293_rep1_RBFOX2','Dan_singleplex_HEK293_rep2_RBFOX2','ENCODE_Dan_singleplex_HEK293_rep1_RBFOX2','ENCODE_Dan_singleplex_HEK293_rep2_RBFOX2']
    
    try:
        os.mkdir('error_files')
    except:
        pass
    R_EXE = config['R_EXE']

    config['GENOME_FA'] = config['GENOMEFA']

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

module region_call:
    snakefile:
        "rules/region_call.py"
    config:
        config

module make_track:
    snakefile:
        config['MAKE_TRACK']
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
        expand("QC/nobarcode_blast_output/{libname}.blast.tsv", libname = libnames),
        expand("output/enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_windows.tsv.gz",
        libname = libnames,
        clip_sample_label = list(set(rbps)-set([config['AS_INPUT']])), # cannot call on itself
        bg_sample_label = [config['AS_INPUT']]
        ,),
        expand("{libname}/bw/{type}/{sample_label}.r1.{strand}.bw", libname = libnames, sample_label = rbps, strand = ['pos', 'neg'], type = ['CITS', 'COV'])

        
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

############# Quality control #################

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

############ region calling ################
# count reads in each region for each library
# line 0 is the {libname}.{sample_label}
use rule partition_bam_reads_r1 from region_call as region_partition_bam_reads with:
    input:
        chrom_size = config['CHROM_SIZES'],
        bam = "{libname}/bams/genome/{sample_label}.genome-mapped.Aligned.sortedByCoord.out.bam",        
        region_partition = config['PARTITION'],
    output:
        counts= "output/counts/genome/vectors/{libname}.{sample_label}.counts",
    params:
        error_out_file = "error_files/{libname}.{sample_label}.partition_bam_reads.err",
        out_file = "stdout/{libname}.{sample_label}.partition_bam_reads.out",
        run_time = "20:00",
        cores = "1",
        memory = "10000",
        job_name = "partition_bam_reads",
        replicate_label = "{libname}.{sample_label}"
    benchmark: "benchmarks/counts/unassigned_experiment.{libname}.{sample_label}.partition_bam_reads.txt"


# concat all the experiments of the same set into table
use rule make_genome_count_table from region_call as make_genome_count_IPonly with:
    input:
        partition=config['PARTITION'],
        replicate_counts = lambda wildcards: expand("output/counts/genome/vectors/{libname}.{sample_label}.counts", 
            libname = [l for l in libnames if wildcards.experiment in l], # TODO: make dictionary
            sample_label = [wildcards.sample_label]),
    output:
        count_table = "output/counts/genome/tables/{experiment}.{sample_label}.tsv.gz",
    params:
        error_out_file = "error_files/{experiment}.{sample_label}.make_count_table.err",
        out_file = "stdout/{experiment}.{sample_label}.make_count_table.out",
        run_time = "00:05:00",
        cores = "1",
        memory = "200",
        job_name = "make_genome_count_table"
    benchmark: "benchmarks/counts/{experiment}.{sample_label}.all_replicates.make_genome_count_table.txt"

use rule calc_partition_nuc from QC

def libname_to_experiment(libname):
    return [e for e in experiments if e in libname][0]

rule fit_clip_betabinom_no_other:
    input:
        nuc = config['PARTITION'].replace(".bed", ".nuc"),
        table = lambda wildcards: "output/counts/genome/tables/"+libname_to_experiment(wildcards.libname)+f".{wildcards.sample_label}.tsv.gz"
    output:
        coef = "output/clip_model_coef/{libname}.{sample_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/clip_distributions/{{experiment}}.{{clip_sample_label}}.{other_label}.clip_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment][wildcards.Input_replicate_label])
    params:
        error_out_file = "error_files/{libname}.{sample_label}.fit_clip_betabinomial_model.err",
        out_file = "stdout/{libname}.{sample_label}.fit_clip_betabinomial_model.out",
        run_time = "1:00:00",
        memory = "1000",
        job_name = "fit_clip_betabinomial_model",
        cores = "1",
    container:
        "docker://algaebrown/beta-binom"
    benchmark: "benchmarks/fit_clip_betabinomial_model/{libname}.{sample_label}.fit_clip.txt"
    shell:
        "{R_EXE} --vanilla {SCRIPT_PATH}/fit_clip_betabinom_no_other_col.R {input.nuc} {input.table} {wildcards.libname} {wildcards.libname}.{wildcards.sample_label} {output.coef}"

rule combine_ip_to_background:
    input:
        count_table = "output/counts/genome/tables/{experiment}.{clip_sample_label}.tsv.gz",
        bg_counts = lambda wildcards: expand("output/counts/genome/vectors/{libname}.{bg_sample_label}.counts", 
            libname = [l for l in libnames if wildcards.experiment in l], # TODO: make dictionary
            bg_sample_label = [wildcards.bg_sample_label])
    output:
        #count_table = "output/counts/genome/tmptables/{experiment}.{clip_sample_label}.tsv",
        combined_count_table = "output/counts/genome/bgtables/{bg_sample_label}/{experiment}.{clip_sample_label}.tsv.gz"
    params:
        error_out_file = "error_files/{experiment}.{bg_sample_label}.{clip_sample_label}.combine.err",
        out_file = "stdout/{experiment}.{bg_sample_label}.{clip_sample_label}.combine.out",
        run_time = "1:00:00",
        cores = "1",
        memory = "1000",
        job_name = "combine_table"
    benchmark: "benchmarks/combine_table/{experiment}.{bg_sample_label}.{clip_sample_label}.combine.txt"
    shell:
        """
        paste <(zcat {input.count_table}) {input.bg_counts}| gzip -c > {output.combined_count_table}
        """

rule call_enriched_windows:
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        accession_rankings = config['ACCESSION_RANKINGS'],
        nuc = config['PARTITION'].replace(".bed", ".nuc"),
        table = lambda wildcards: f"output/counts/genome/bgtables/{wildcards.bg_sample_label}/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
        parameters = lambda wildcards: "output/clip_model_coef/{libname}.{bg_sample_label}.tsv",
    output:
        "output/threshold_scan/{libname}.{clip_sample_label}.{bg_sample_label}.threshold_data.tsv",
        "output/tested_windows/{libname}.{clip_sample_label}.{bg_sample_label}.tested_windows.tsv.gz",
        "output/enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_windows.tsv.gz",
        "output/enrichment_summaries/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_window_feature_data.tsv",
        "output/enrichment_summaries/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_window_transcript_data.tsv",
        "output/enrichment_summaries/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_window_gene_data.tsv",
        "output/all_reads/{libname}.{clip_sample_label}.{bg_sample_label}.all_reads_fractions_feature_data.tsv",
        "output/all_reads/{libname}.{clip_sample_label}.{bg_sample_label}.all_reads_odds_feature_data.tsv",
        "output/all_reads/{libname}.{clip_sample_label}.{bg_sample_label}.all_reads_odds_transcript_data.tsv",
        "output/all_reads/{libname}.{clip_sample_label}.{bg_sample_label}.all_reads_odds_feature_gc_data.tsv",
        "output/figures/threshold_scan/{libname}.{clip_sample_label}.{bg_sample_label}.threshold_scan.pdf",
        "output/figures/enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_window_coverage.pdf",
        "output/figures/enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_window_rates.pdf",
        "output/figures/enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_window_counts.linear.pdf",
        "output/figures/enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_window_counts.log10.pdf",
        "output/figures/enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_window_odds.feature.pdf",
        "output/figures/enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_window_odds.all_transcript_types.pdf",
        "output/figures/enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_window_odds.select_transcript_types.pdf",
        "output/figures/enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.enriched_window_counts.per_gene_feature.pdf",
        "output/figures/all_reads/{libname}.{clip_sample_label}.{bg_sample_label}.all_reads_fractions.feature.pdf",
        "output/figures/all_reads/{libname}.{clip_sample_label}.{bg_sample_label}.all_reads_odds.feature.pdf",
        "output/figures/all_reads/{libname}.{clip_sample_label}.{bg_sample_label}.all_reads_odds.all_transcript_types.pdf",
        "output/figures/all_reads/{libname}.{clip_sample_label}.{bg_sample_label}.all_reads_odds.feature_gc.pdf"
    params:
        error_out_file = "error_files/{libname}.{clip_sample_label}.{bg_sample_label}.call_enriched_windows.err",
        out_file = "stdout/{libname}.{clip_sample_label}.{bg_sample_label}.call_enriched_windows.out",
        run_time = "00:25:00",
        memory = "1000",
        job_name = "call_enriched_windows",
        cores = "1",
    benchmark: "benchmarks/call_enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.call_enriched_windows.txt"
    # container:
    #     "docker://algaebrown/beta-binom" # TODO: THIS FUCKING SHIT WORKS WITH COPY AND PASTE BUT NOT SNAKEMAKE. no error msg
    shell:
        """
        {R_EXE} --vanilla {SCRIPT_PATH}/call_enriched_windows.R \
            {input.nuc} \
            {input.table} \
            {input.accession_rankings} \
            {input.feature_annotations} \
            {input.parameters} \
            {wildcards.libname}.{wildcards.bg_sample_label} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {wildcards.libname}.{wildcards.clip_sample_label}.{wildcards.bg_sample_label}
        """


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
        cores = 1,
    shell:
        """
        module load ucsctools
        bedGraphToBigWig {input.bedgraph} {params.chr_size} {output.bw}
        """

