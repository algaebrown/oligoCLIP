import pandas as pd
R_EXE = config['R_EXE']
BLACKLIST = config['BLACKLIST'] if 'BLACKLIST' in config else None
rbps = config['rbps']
experiments = config['experiments']
libnames = config['libnames']
SCRIPT_PATH=config['SCRIPT_PATH']
manifest = pd.read_table(config['MANIFEST'], index_col = False, sep = ',')
UNINFORMATIVE_READ = 3 - int(config['INFORMATIVE_READ']) # whether read 1 or read 2 is informative

try:
    os.mkdir('log')
except:
    pass

module region_call:
    snakefile:
        "region_call.py"
    config:
        config

module QC:
    snakefile:
        "QC.py"
    config:
        config

def libname_to_experiment(libname):
    return manifest.loc[manifest['libname']==libname, 'experiment'].iloc[0]
def experiment_to_libname(experiment):
    libnames = manifest.loc[manifest['experiment']==experiment, 'libname'].tolist()
    assert len(libnames)>0
    return libnames
############ region calling ################
# count reads in each region for each library
# line 0 is the {libname}.{sample_label}
use rule partition_bam_reads from region_call as region_partition_bam_reads with:
    input:
        chrom_size = config['CHROM_SIZES'],
        bam = "{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",        
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
            libname = experiment_to_libname(wildcards.experiment), # TODO: make dictionary
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

################# for IgG normalization only ##########################
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
        log = "log/{libname}.{sample_label}.fit_binom.log"
    container:
        "docker://algaebrown/beta-binom"
    benchmark: "benchmarks/fit_clip_betabinomial_model/{libname}.{sample_label}.fit_clip.txt"
    shell:
        """
        {R_EXE} --vanilla {SCRIPT_PATH}fit_clip_betabinom_no_other_col.R {input.nuc} {input.table} {wildcards.libname} {wildcards.libname}.{wildcards.sample_label} {output.coef} &> {params.log}
        """

rule combine_ip_to_background:
    input:
        count_table = "output/counts/genome/tables/{experiment}.{clip_sample_label}.tsv.gz",
        bg_counts = lambda wildcards: expand("output/counts/genome/vectors/{libname}.{bg_sample_label}.counts", 
            libname = experiment_to_libname(wildcards.experiment), # TODO: make dictionary
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
            {wildcards.libname}.{wildcards.clip_sample_label}.{wildcards.bg_sample_label} &> {params.out_file}
        """

####################### internal normalization ######################
rule sum_all_other_background:
    input:
        lambda wildcards: expand("output/counts/genome/vectors/{libname}.{sample_label}.counts",
        libname = ["{libname}"],
        sample_label = list(set(rbps)-set([wildcards.clip_sample_label])-set([config['AS_INPUT']]))
        )
    output:
        counts= "internal_output/counts/genome/vectors/{libname}.{clip_sample_label}.counts",
    params:
        error_out_file = "error_files/{libname}.{clip_sample_label}.sum_read.err",
        out_file = "stdout/{libname}.{clip_sample_label}.sum_read.out",
        run_time = "20:00",
        cores = "1",
        memory = "10000",
        job_name = "sum_other_libs",
        replicate_label = "{libname}.internal"
    benchmark: "benchmarks/counts/unassigned_experiment.{libname}.{clip_sample_label}.sum_read.txt"
    shell:
        """
        awk '{{arr[FNR]+=$1}}END{{for(i=2;i<=FNR;i+=1){{print arr[i]}} }}' {input} | awk 'BEGIN {{print \"{params.replicate_label}\"}} {{print}}' > {output}
        """

rule make_internal_genome_count_table:
    input:
        partition=config['PARTITION'],
        replicate_counts = lambda wildcards: expand("internal_output/counts/genome/vectors/{libname}.{clip_sample_label}.counts", 
            libname = experiment_to_libname(wildcards.experiment), # TODO: make dictionary
            clip_sample_label = [wildcards.clip_sample_label]),
    output:
        count_table = "internal_output/counts/genome/tables/{experiment}.{clip_sample_label}.tsv.gz",
    params:
        error_out_file = "error_files/{experiment}.{clip_sample_label}.bgsum.make_count_table.err",
        out_file = "stdout/{experiment}.{clip_sample_label}.bgsum.make_count_table.out",
        run_time = "00:05:00",
        cores = "1",
        memory = "200",
        job_name = "make_genome_count_table"
    benchmark: "benchmarks/counts/{experiment}.{clip_sample_label}.bgsum.all_replicates.make_genome_count_table.txt"
    shell:
        "paste <(awk -v OFS=\"\\t\" 'BEGIN {{print \"chr\\tstart\\tend\\tname\\tscore\\tstrand\"}} {{print}}' {input.partition}) {input.replicate_counts} | gzip -c > {output.count_table};"

rule fit_clip_betabinom_internal:
    # overdispersion from that sum of other CLIP library
    input:
        nuc = config['PARTITION'].replace(".bed", ".nuc"),
        table = lambda wildcards: "internal_output/counts/genome/tables/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz"
    output:
        coef = "internal_output/clip_model_coef/{libname}.{clip_sample_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/clip_distributions/{{experiment}}.{{clip_sample_label}}.{other_label}.clip_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment][wildcards.Input_replicate_label])
    params:
        error_out_file = "error_files/{libname}.{clip_sample_label}.fit_internal_clip_betabinomial_model.err",
        out_file = "stdout/{libname}.{clip_sample_label}.fit_internal_clip_betabinomial_model.out",
        run_time = "1:00:00",
        memory = "10000",
        job_name = "fit_internal_clip_betabinomial_model",
        cores = "1",
    container:
        "docker://algaebrown/beta-binom"
    benchmark: "benchmarks/fit_clip_betabinomial_model/{libname}.{clip_sample_label}.fit_clip.txt"
    shell:
        """
        {R_EXE} --vanilla {SCRIPT_PATH}/fit_clip_betabinom_no_other_col.R {input.nuc} {input.table} {wildcards.libname} {wildcards.libname}.internal {output.coef} 2> {params.out_file}
        """

rule combine_ip_to_internal_background: # bg_sample_label = 'internal'
    input:
        count_table = "output/counts/genome/tables/{experiment}.{clip_sample_label}.tsv.gz",
        bg_counts = lambda wildcards: expand("internal_output/counts/genome/vectors/{libname}.{clip_sample_label}.counts", 
            libname = experiment_to_libname(wildcards.experiment), # TODO: make dictionary
            clip_sample_label = [wildcards.clip_sample_label])
    output:
        #count_table = "output/counts/genome/tmptables/{experiment}.{clip_sample_label}.tsv",
        combined_count_table = "internal_output/counts/genome/bgtables/internal/{experiment}.{clip_sample_label}.tsv.gz"
    params:
        error_out_file = "error_files/{experiment}.internal.{clip_sample_label}.combine.err",
        out_file = "stdout/{experiment}.internal.{clip_sample_label}.combine.out",
        run_time = "1:00:00",
        cores = "1",
        memory = "10000",
        job_name = "combine_table"
    benchmark: "benchmarks/combine_table/{experiment}.internal.{clip_sample_label}.combine.txt"
    shell:
        """
        paste <(zcat {input.count_table}) {input.bg_counts}| gzip -c > {output.combined_count_table}
        """
rule call_enriched_windows_internal:
    # overdispersion from parameters
    # mean from table, sum of all other CLIP libraries
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        accession_rankings = config['ACCESSION_RANKINGS'],
        nuc = config['PARTITION'].replace(".bed", ".nuc"),
        table = lambda wildcards: f"internal_output/counts/genome/bgtables/internal/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
        #parameters = lambda wildcards: "internal_output/clip_model_coef/{libname}.{clip_sample_label}.tsv", # overdispersion from that sum of other CLIP library
        parameters = lambda wildcards: "output/clip_model_coef/{libname}.{clip_sample_label}.tsv", # overdispersion from that CLIP library
    output:
        "internal_output/threshold_scan/{libname}.{clip_sample_label}.threshold_data.tsv",
        "internal_output/tested_windows/{libname}.{clip_sample_label}.tested_windows.tsv.gz",
        "internal_output/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
        "internal_output/enrichment_summaries/{libname}.{clip_sample_label}.enriched_window_feature_data.tsv",
        "internal_output/enrichment_summaries/{libname}.{clip_sample_label}.enriched_window_transcript_data.tsv",
        "internal_output/enrichment_summaries/{libname}.{clip_sample_label}.enriched_window_gene_data.tsv",
        "internal_output/all_reads/{libname}.{clip_sample_label}.all_reads_fractions_feature_data.tsv",
        "internal_output/all_reads/{libname}.{clip_sample_label}.all_reads_odds_feature_data.tsv",
        "internal_output/all_reads/{libname}.{clip_sample_label}.all_reads_odds_transcript_data.tsv",
        "internal_output/all_reads/{libname}.{clip_sample_label}.all_reads_odds_feature_gc_data.tsv",
        "internal_output/figures/threshold_scan/{libname}.{clip_sample_label}.threshold_scan.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_coverage.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_rates.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_counts.linear.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_counts.log10.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_odds.feature.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_odds.all_transcript_types.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_odds.select_transcript_types.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_counts.per_gene_feature.pdf",
        "internal_output/figures/all_reads/{libname}.{clip_sample_label}.all_reads_fractions.feature.pdf",
        "internal_output/figures/all_reads/{libname}.{clip_sample_label}.all_reads_odds.feature.pdf",
        "internal_output/figures/all_reads/{libname}.{clip_sample_label}.all_reads_odds.all_transcript_types.pdf",
        "internal_output/figures/all_reads/{libname}.{clip_sample_label}.all_reads_odds.feature_gc.pdf"
    params:
        error_out_file = "error_files/{libname}.{clip_sample_label}.internal.call_enriched_windows.err",
        out_file = "stdout/{libname}.{clip_sample_label}.internal.call_enriched_windows.out",
        run_time = "00:25:00",
        memory = "10000",
        job_name = "call_enriched_windows",
        cores = "8",
    benchmark: "benchmarks/call_enriched_windows/{libname}.{clip_sample_label}.internal.call_enriched_windows.txt"
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
            {wildcards.libname}.internal \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {wildcards.libname}.{wildcards.clip_sample_label} &> {params.out_file}
        """

rule call_enriched_windows_internal_dlogs:
    # overdispersion from parameters
    # mean from table, sum of all other CLIP libraries
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        accession_rankings = config['ACCESSION_RANKINGS'],
        nuc = config['PARTITION'].replace(".bed", ".nuc"),
        table = lambda wildcards: f"internal_output/counts/genome/bgtables/internal/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
        #parameters = lambda wildcards: "internal_output/clip_model_coef/{libname}.{clip_sample_label}.tsv", # overdispersion from that sum of other CLIP library
        parameters = lambda wildcards: "output/clip_model_coef/{libname}.{clip_sample_label}.tsv", # overdispersion from that CLIP library
    output:
        "internal_output/threshold_scan/{libname}.{clip_sample_label}.dlogs.threshold_data.tsv",
        "internal_output/tested_windows/{libname}.{clip_sample_label}.dlogs.tested_windows.tsv.gz",
        "internal_output/enriched_windows/{libname}.{clip_sample_label}.dlogs.enriched_windows.tsv.gz",
        "internal_output/enrichment_summaries/{libname}.{clip_sample_label}.dlogs.enriched_window_feature_data.tsv",
        "internal_output/enrichment_summaries/{libname}.{clip_sample_label}.dlogs.enriched_window_transcript_data.tsv",
        "internal_output/enrichment_summaries/{libname}.{clip_sample_label}.dlogs.enriched_window_gene_data.tsv",
        "internal_output/all_reads/{libname}.{clip_sample_label}.dlogs.all_reads_fractions_feature_data.tsv",
        "internal_output/all_reads/{libname}.{clip_sample_label}.dlogs.all_reads_odds_feature_data.tsv",
        "internal_output/all_reads/{libname}.{clip_sample_label}.dlogs.all_reads_odds_transcript_data.tsv",
        "internal_output/all_reads/{libname}.{clip_sample_label}.dlogs.all_reads_odds_feature_gc_data.tsv",
        "internal_output/figures/threshold_scan/{libname}.{clip_sample_label}.dlogs.threshold_scan.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.dlogs.enriched_window_coverage.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.dlogs.enriched_window_rates.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.dlogs.enriched_window_counts.linear.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.dlogs.enriched_window_counts.log10.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.dlogs.enriched_window_odds.feature.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.dlogs.enriched_window_odds.all_transcript_types.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.dlogs.enriched_window_odds.select_transcript_types.pdf",
        "internal_output/figures/enriched_windows/{libname}.{clip_sample_label}.dlogs.enriched_window_counts.per_gene_feature.pdf",
        "internal_output/figures/all_reads/{libname}.{clip_sample_label}.dlogs.all_reads_fractions.feature.pdf",
        "internal_output/figures/all_reads/{libname}.{clip_sample_label}.dlogs.all_reads_odds.feature.pdf",
        "internal_output/figures/all_reads/{libname}.{clip_sample_label}.dlogs.all_reads_odds.all_transcript_types.pdf",
        "internal_output/figures/all_reads/{libname}.{clip_sample_label}.dlogs.all_reads_odds.feature_gc.pdf"
    params:
        error_out_file = "error_files/{libname}.{clip_sample_label}.dlogs.internal.call_enriched_windows.err",
        out_file = "stdout/{libname}.{clip_sample_label}.dlogs.internal.call_enriched_windows.out",
        run_time = "00:25:00",
        memory = "10000",
        job_name = "call_enriched_windows",
        cores = "8",
    benchmark: "benchmarks/call_enriched_windows/{libname}.{clip_sample_label}.dlogs.internal.call_enriched_windows.txt"
    # container:
    #     "docker://algaebrown/beta-binom" # TODO: THIS FUCKING SHIT WORKS WITH COPY AND PASTE BUT NOT SNAKEMAKE. no error msg
    shell:
        """
        {R_EXE} --vanilla {SCRIPT_PATH}/call_enriched_windows_threshold_d_log.R \
            {input.nuc} \
            {input.table} \
            {input.accession_rankings} \
            {input.feature_annotations} \
            {input.parameters} \
            {wildcards.libname}.internal \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {wildcards.libname}.{wildcards.clip_sample_label}.dlogs &> {params.out_file}
        """


rule find_reproducible_enriched_windows:
    input: # "internal_output/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
        windows = lambda wildcards: expand(
            "internal_output/enriched_windows/{libname}.{{clip_sample_label}}.enriched_windows.tsv.gz", 
            libname = manifest.loc[manifest['experiment']==wildcards.experiment, 'libname'].tolist()) # find all library of the same experiment
    output:
        reproducible_windows = "output/reproducible_enriched_windows/{experiment}.{clip_sample_label}.reproducible_enriched_windows.tsv.gz",
        linear_bar = "output/figures/reproducible_enriched_windows/{experiment}.{clip_sample_label}.reproducible_enriched_window_counts.linear.pdf",
        log_bar = "output/figures/reproducible_enriched_windows/{experiment}.{clip_sample_label}.reproducible_enriched_window_counts.log10.pdf"
    params:
        error_file = "stderr/{experiment}.{clip_sample_label}.find_reproducible_enriched_windows.err",
        out_file = "stdout/{experiment}.{clip_sample_label}.find_reproducible_enriched_windows.out",
        run_time = "5:00",
        memory = "2000",
        job_name = "find_reproducible_enriched_windows"
    benchmark: "benchmarks/find_reproducible_enriched_windows/{experiment}.{clip_sample_label}.all_replicates.reproducible.txt"
    shell:
        "{R_EXE} --vanilla {SCRIPT_PATH}/identify_reproducible_windows.R internal_output/enriched_windows/ {wildcards.clip_sample_label} " + (BLACKLIST if BLACKLIST is not None else "") 
