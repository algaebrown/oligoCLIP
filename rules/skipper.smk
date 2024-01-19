import pandas as pd
locals().update(config)
rbps = config['rbps']
experiments = config['experiments']
libnames = config['libnames']
manifest = pd.read_table(config['MANIFEST'], index_col = False, sep = ',')


def libname_to_experiment(libname):
    return manifest.loc[manifest['libname']==libname, 'experiment'].iloc[0]
def experiment_to_libname(experiment):
    libnames = manifest.loc[manifest['experiment']==experiment, 'libname'].tolist()
    assert len(libnames)>0
    return libnames
############ counting ################
# count reads in each region for each library
# line 0 is the {libname}.{sample_label}
rule partition_bam_reads:
    input:
        chrom_size = config['CHROM_SIZES'],
        bam = "{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",        
        region_partition = config['PARTITION'],
    output:
        counts= "counts/genome/vectors/{libname}.{sample_label}.counts",
    params:
        error_out_file = "error_files/partition_bam_reads.{libname}.{sample_label}.err",
        out_file = "stdout/partition_bam_reads.{libname}.{sample_label}.out",
        run_time = "20:00",
        cores = "1",
        memory = 10000,
        replicate_label = "{libname}.{sample_label}",
        uninformative_read = config['UNINFORMATIVE_READ']
    benchmark: "benchmarks/counts/unassigned_experiment.{libname}.{sample_label}.partition_bam_reads.txt"
    conda:
        "envs/bedtools.yaml"
    shell:
        "bedtools bamtobed -i {input.bam} | awk '($1 != \"chrEBV\") && ($4 !~ \"/{params.uninformative_read}$\")' | bedtools flank -s -l 1 -r 0 -g {input.chrom_size} -i - | bedtools shift -p -1 -m 1 -g {input.chrom_size} -i - | bedtools coverage -counts -s -a {input.region_partition} -b - | cut -f 7 | awk 'BEGIN {{print \"{params.replicate_label}\"}} {{print}}' > {output.counts};"

# concat all reps of the same experiment into 1 table with annotation
# outputs columns: [annotation] [repcounts]
rule make_genome_count_table:
    input:
        partition=config['PARTITION'],
        replicate_counts = lambda wildcards: expand("counts/genome/vectors/{libname}.{sample_label}.counts", 
            libname = experiment_to_libname(wildcards.experiment), # TODO: make dictionary
            sample_label = [wildcards.sample_label]),
    output:
        count_table = "counts/genome/tables/{experiment}.{sample_label}.tsv.gz",
    params:
        error_out_file = "error_files/{experiment}.{sample_label}.make_count_table.err",
        out_file = "stdout/{experiment}.{sample_label}.make_count_table.out",
        run_time = "00:05:00",
        cores = "1",
        memory = 200,
    benchmark: "benchmarks/counts/{experiment}.{sample_label}.all_replicates.make_genome_count_table.txt"
    shell:
        "paste <(zcat {input.partition} | awk -v OFS=\"\\t\" 'BEGIN {{print \"chr\\tstart\\tend\\tname\\tscore\\tstrand\"}} {{print $1,$2,$3,$4,$5,$6}}' ) {input.replicate_counts} | gzip -c > {output.count_table}"

rule calc_partition_nuc:
    input:
        partition = config['PARTITION'],
        genome = config['GENOMEFA']
    output:
        nuc = config['PARTITION'].replace(".bed", ".nuc")
    params:
        error_out_file = "stderr/calc_partition_nuc.err",
        out_file = "stdout/calc_partition_nuc.out",
        run_time = "00:10:00",
        memory = 1000,
    benchmark: "benchmarks/partition_nuc.txt"
    container:
        "docker://howardxu520/skipper:bigwig_1.0"
    shell:
        "bedtools nuc -s -fi {input.genome} -bed {input.partition} > {output.nuc}"

rule fit_clip_betabinom:
    input:
        nuc = ancient(config['PARTITION'].replace(".bed", ".nuc")),
        table = lambda wildcards: "counts/genome/tables/"+libname_to_experiment(wildcards.libname)+f".{wildcards.sample_label}.tsv.gz"
    output:
        coef = "skipper/clip_model_coef/{libname}.{sample_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/clip_distributions/{{experiment}}.{{clip_sample_label}}.{other_label}.clip_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment][wildcards.Input_replicate_label])
    params:
        error_out_file = "error_files/fit_clip_betabinomial_model.{libname}.{sample_label}.err",
        out_file = "stdout/fit_clip_betabinomial_model.{libname}.{sample_label}.out",
        run_time = "1:00:00",
        memory = 40000,
        cores = "1"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    benchmark: "benchmarks/fit_clip_betabinomial_model/{libname}.{sample_label}.fit_clip.txt"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}fit_clip_betabinom_no_other_col.R {input.nuc} {input.table} {wildcards.libname} {wildcards.libname}.{wildcards.sample_label} {output.coef} 
        """

rule combine_ip_to_background:
    input:
        count_table = "counts/genome/tables/{experiment}.{clip_sample_label}.tsv.gz",
        bg_counts = lambda wildcards: expand("counts/genome/vectors/{libname}.{bg_sample_label}.counts", 
            libname = experiment_to_libname(wildcards.experiment), # TODO: make dictionary
            bg_sample_label = [wildcards.bg_sample_label])
    output:
        combined_count_table = "counts/genome/bgtables/{bg_sample_label}/{experiment}.{clip_sample_label}.tsv.gz"
    params:
        error_out_file = "error_files/combine.{experiment}.{bg_sample_label}.{clip_sample_label}",
        out_file = "stdout/combine.{experiment}.{bg_sample_label}.{clip_sample_label}",
        run_time = "1:00:00",
        cores = "1",
        memory = 1000,
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
        table = lambda wildcards: f"counts/genome/bgtables/{wildcards.bg_sample_label}/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
        parameters = lambda wildcards: "skipper/clip_model_coef/{libname}.{bg_sample_label}.tsv",
    output:
        "skipper/{bg_sample_label}/threshold_scan/{libname}.{clip_sample_label}.threshold_data.tsv",
        "skipper/{bg_sample_label}/tested_windows/{libname}.{clip_sample_label}.tested_windows.tsv.gz",
        "skipper/{bg_sample_label}/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
        "skipper/{bg_sample_label}/enrichment_summaries/{libname}.{clip_sample_label}.enriched_window_feature_data.tsv",
        "skipper/{bg_sample_label}/enrichment_summaries/{libname}.{clip_sample_label}.enriched_window_transcript_data.tsv",
        "skipper/{bg_sample_label}/enrichment_summaries/{libname}.{clip_sample_label}.enriched_window_gene_data.tsv",
        "skipper/{bg_sample_label}/all_reads/{libname}.{clip_sample_label}.all_reads_fractions_feature_data.tsv",
        "skipper/{bg_sample_label}/all_reads/{libname}.{clip_sample_label}.all_reads_odds_feature_data.tsv",
        "skipper/{bg_sample_label}/all_reads/{libname}.{clip_sample_label}.all_reads_odds_transcript_data.tsv",
        "skipper/{bg_sample_label}/all_reads/{libname}.{clip_sample_label}.all_reads_odds_feature_gc_data.tsv",
        "skipper/{bg_sample_label}/figures/threshold_scan/{libname}.{clip_sample_label}.threshold_scan.pdf",
        "skipper/{bg_sample_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_coverage.pdf",
        "skipper/{bg_sample_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_rates.pdf",
        "skipper/{bg_sample_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_counts.linear.pdf",
        "skipper/{bg_sample_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_counts.log10.pdf",
        "skipper/{bg_sample_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_odds.feature.pdf",
        "skipper/{bg_sample_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_odds.all_transcript_types.pdf",
        "skipper/{bg_sample_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_odds.select_transcript_types.pdf",
        "skipper/{bg_sample_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_counts.per_gene_feature.pdf",
        "skipper/{bg_sample_label}/figures/all_reads/{libname}.{clip_sample_label}.all_reads_fractions.feature.pdf",
        "skipper/{bg_sample_label}/figures/all_reads/{libname}.{clip_sample_label}.all_reads_odds.feature.pdf",
        "skipper/{bg_sample_label}/figures/all_reads/{libname}.{clip_sample_label}.all_reads_odds.all_transcript_types.pdf",
        "skipper/{bg_sample_label}/figures/all_reads/{libname}.{clip_sample_label}.all_reads_odds.feature_gc.pdf"
    params:
        error_out_file = "error_files/call_enriched_window.bg.{libname}.{clip_sample_label}.{bg_sample_label}.err",
        out_file = "stdout/call_enriched_window.bg.{libname}.{clip_sample_label}.{bg_sample_label}.out",
        run_time = "00:25:00",
        memory = 80000,
        cores = "1",
        root_folder = lambda wildcards, output: str(output[0]).split('threshold_scan')[0]
    benchmark: "benchmarks/call_enriched_windows/{libname}.{clip_sample_label}.{bg_sample_label}.call_enriched_windows.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/call_enriched_windows.R \
            {input.nuc} \
            {input.table} \
            {input.accession_rankings} \
            {input.feature_annotations} \
            {input.parameters} \
            {wildcards.libname}.{wildcards.bg_sample_label} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {wildcards.libname}.{wildcards.clip_sample_label}.{wildcards.bg_sample_label} \
            {params.root_folder}
        """

####################### CC: complementary control ######################
rule sum_all_other_background_as_CC:
    input:
        lambda wildcards: expand("counts/genome/vectors/{libname}.{sample_label}.counts",
        libname = ["{libname}"],
        sample_label = list(set(rbps)-set([wildcards.clip_sample_label])-set(config['AS_INPUT']))
        )
    output:
        counts= "counts_CC/genome/vectors/{libname}.{clip_sample_label}.counts",
    params:
        error_out_file = "error_files/sum_reads.{libname}.{clip_sample_label}",
        out_file = "stdout/sum_reads.{libname}.{clip_sample_label}",
        run_time = "20:00",
        cores = "1",
        memory = 10000,
        replicate_label = "{libname}.internal"
    benchmark: "benchmarks/counts/unassigned_experiment.{libname}.{clip_sample_label}.sum_read.txt"
    shell:
        """
        awk '{{arr[FNR]+=$1}}END{{for(i=2;i<=FNR;i+=1){{print arr[i]}} }}' {input} | awk 'BEGIN {{print \"{params.replicate_label}\"}} {{print}}' > {output}
        """

rule combine_ip_to_CC: # bg_sample_label = 'internal'
    input:
        count_table = "counts/genome/tables/{experiment}.{clip_sample_label}.tsv.gz",
        bg_counts = lambda wildcards: expand("counts_CC/genome/vectors/{libname}.{clip_sample_label}.counts", 
            libname = experiment_to_libname(wildcards.experiment), # TODO: make dictionary
            clip_sample_label = [wildcards.clip_sample_label])
    output:
        combined_count_table = "counts_CC/genome/bgtables/internal/{experiment}.{clip_sample_label}.tsv.gz"
    params:
        error_out_file = "error_files/combine.CC.{experiment}.{clip_sample_label}.err",
        out_file = "stdout/combine.CC{experiment}.{clip_sample_label}.out",
        run_time = "1:00:00",
        cores = "1",
        memory = 10000,
    benchmark: "benchmarks/combine_table/{experiment}.internal.{clip_sample_label}.combine.txt"
    shell:
        """
        paste <(zcat {input.count_table}) {input.bg_counts}| gzip -c > {output.combined_count_table}
        """
rule call_enriched_windows_CC:
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        accession_rankings = config['ACCESSION_RANKINGS'],
        nuc = config['PARTITION'].replace(".bed", ".nuc"),
        table = lambda wildcards: f"counts_CC/genome/bgtables/internal/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
        parameters = lambda wildcards: "skipper/clip_model_coef/{libname}.{clip_sample_label}.tsv", # overdispersion from that CLIP library
    output:
        "skipper_CC/threshold_scan/{libname}.{clip_sample_label}.threshold_data.tsv",
        "skipper_CC/tested_windows/{libname}.{clip_sample_label}.tested_windows.tsv.gz",
        "skipper_CC/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
        "skipper_CC/enrichment_summaries/{libname}.{clip_sample_label}.enriched_window_feature_data.tsv",
        "skipper_CC/enrichment_summaries/{libname}.{clip_sample_label}.enriched_window_transcript_data.tsv",
        "skipper_CC/enrichment_summaries/{libname}.{clip_sample_label}.enriched_window_gene_data.tsv",
        "skipper_CC/all_reads/{libname}.{clip_sample_label}.all_reads_fractions_feature_data.tsv",
        "skipper_CC/all_reads/{libname}.{clip_sample_label}.all_reads_odds_feature_data.tsv",
        "skipper_CC/all_reads/{libname}.{clip_sample_label}.all_reads_odds_transcript_data.tsv",
        "skipper_CC/all_reads/{libname}.{clip_sample_label}.all_reads_odds_feature_gc_data.tsv",
        "skipper_CC/figures/threshold_scan/{libname}.{clip_sample_label}.threshold_scan.pdf",
        "skipper_CC/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_coverage.pdf",
        "skipper_CC/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_rates.pdf",
        "skipper_CC/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_counts.linear.pdf",
        "skipper_CC/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_counts.log10.pdf",
        "skipper_CC/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_odds.feature.pdf",
        "skipper_CC/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_odds.all_transcript_types.pdf",
        "skipper_CC/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_odds.select_transcript_types.pdf",
        "skipper_CC/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_counts.per_gene_feature.pdf",
        "skipper_CC/figures/all_reads/{libname}.{clip_sample_label}.all_reads_fractions.feature.pdf",
        "skipper_CC/figures/all_reads/{libname}.{clip_sample_label}.all_reads_odds.feature.pdf",
        "skipper_CC/figures/all_reads/{libname}.{clip_sample_label}.all_reads_odds.all_transcript_types.pdf",
        "skipper_CC/figures/all_reads/{libname}.{clip_sample_label}.all_reads_odds.feature_gc.pdf"
    params:
        error_out_file = "error_files/call_enriched_windows.CC.{libname}.{clip_sample_label}",
        out_file = "stdout/call_enriched_windows.CC.{libname}.{clip_sample_label}",
        run_time = "00:25:00",
        memory = 40000,
        cores = "8",
        root_folder = 'skipper_CC'
    benchmark: "benchmarks/call_enriched_windows/{libname}.{clip_sample_label}.internal.call_enriched_windows.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/call_enriched_windows.R \
            {input.nuc} \
            {input.table} \
            {input.accession_rankings} \
            {input.feature_annotations} \
            {input.parameters} \
            {wildcards.libname}.internal \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {params.root_folder}

        """

####################### external normalization ######################
use rule partition_bam_reads as partition_external_bams with:
    input:
        chrom_size = config['CHROM_SIZES'],
        bam = lambda wildcards: ancient(config['external_bam'][wildcards.external_label]['file']),   
        region_partition = config['PARTITION'],
    output:
        counts= "counts/genome/vectors/external.{external_label}.counts",
    params:
        error_out_file = "error_files/partition_bam_reads.external.{external_label}.err",
        out_file = "stdout/partition_bam_reads.external.{external_label}.out",
        run_time = "20:00",
        cores = "1",
        memory = 10000,
        replicate_label = "external.{external_label}",
        uninformative_read = lambda wildcards: 3 - config['external_bam'][wildcards.external_label]['INFORMATIVE_READ']
    benchmark: "benchmarks/counts/unassigned_experiment.external.{external_label}.partition_bam_reads.txt"

use rule combine_ip_to_CC as combine_ip_to_external with:
    input:
        count_table = "counts/genome/tables/{experiment}.{clip_sample_label}.tsv.gz",
        bg_counts = rules.partition_external_bams.output.counts
    output:
        combined_count_table = "counts_external/genome/{external_label}/{experiment}.{clip_sample_label}.tsv.gz"
    params:
        error_out_file = "error_files/combine.{external_label}.{experiment}.{clip_sample_label}.err",
        out_file = "stdout/combine.{external_label}.{experiment}.{clip_sample_label}.out",
        run_time = "1:00:00",
        cores = "1",
        memory = "10000",
    benchmark: "benchmarks/combine_table/{experiment}.{external_label}.{clip_sample_label}.combine.txt"

rule call_enriched_windows_external:
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        accession_rankings = config['ACCESSION_RANKINGS'],
        nuc = config['PARTITION'].replace(".bed", ".nuc"),
        table = lambda wildcards: f"counts_external/genome/{wildcards.external_label}/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
        parameters = lambda wildcards: "skipper/clip_model_coef/{libname}.{clip_sample_label}.tsv", # overdispersion from that CLIP library
    output:
        "skipper_external/{external_label}/threshold_scan/{libname}.{clip_sample_label}.threshold_data.tsv",
        "skipper_external/{external_label}/tested_windows/{libname}.{clip_sample_label}.tested_windows.tsv.gz",
        "skipper_external/{external_label}/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
        "skipper_external/{external_label}/enrichment_summaries/{libname}.{clip_sample_label}.enriched_window_feature_data.tsv",
        "skipper_external/{external_label}/enrichment_summaries/{libname}.{clip_sample_label}.enriched_window_transcript_data.tsv",
        "skipper_external/{external_label}/enrichment_summaries/{libname}.{clip_sample_label}.enriched_window_gene_data.tsv",
        "skipper_external/{external_label}/all_reads/{libname}.{clip_sample_label}.all_reads_fractions_feature_data.tsv",
        "skipper_external/{external_label}/all_reads/{libname}.{clip_sample_label}.all_reads_odds_feature_data.tsv",
        "skipper_external/{external_label}/all_reads/{libname}.{clip_sample_label}.all_reads_odds_transcript_data.tsv",
        "skipper_external/{external_label}/all_reads/{libname}.{clip_sample_label}.all_reads_odds_feature_gc_data.tsv",
        "skipper_external/{external_label}/figures/threshold_scan/{libname}.{clip_sample_label}.threshold_scan.pdf",
        "skipper_external/{external_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_coverage.pdf",
        "skipper_external/{external_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_rates.pdf",
        "skipper_external/{external_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_counts.linear.pdf",
        "skipper_external/{external_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_counts.log10.pdf",
        "skipper_external/{external_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_odds.feature.pdf",
        "skipper_external/{external_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_odds.all_transcript_types.pdf",
        "skipper_external/{external_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_odds.select_transcript_types.pdf",
        "skipper_external/{external_label}/figures/enriched_windows/{libname}.{clip_sample_label}.enriched_window_counts.per_gene_feature.pdf",
        "skipper_external/{external_label}/figures/all_reads/{libname}.{clip_sample_label}.all_reads_fractions.feature.pdf",
        "skipper_external/{external_label}/figures/all_reads/{libname}.{clip_sample_label}.all_reads_odds.feature.pdf",
        "skipper_external/{external_label}/figures/all_reads/{libname}.{clip_sample_label}.all_reads_odds.all_transcript_types.pdf",
        "skipper_external/{external_label}/figures/all_reads/{libname}.{clip_sample_label}.all_reads_odds.feature_gc.pdf"
    params:
        error_out_file = "error_files/call_enriched_windows.{external_label}.{libname}.{clip_sample_label}",
        out_file = "stdout/call_enriched_windows.{external_label}.{libname}.{clip_sample_label}",
        run_time = "00:25:00",
        memory = 40000,
        cores = "8",
        root_folder = 'skipper_external/{external_label}'
    benchmark: "benchmarks/call_enriched_windows/{libname}.{clip_sample_label}.{external_label}.call_enriched_windows.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/call_enriched_windows.R \
            {input.nuc} \
            {input.table} \
            {input.accession_rankings} \
            {input.feature_annotations} \
            {input.parameters} \
            external.{wildcards.external_label} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {params.root_folder}

        """


# rule find_reproducible_enriched_windows:
#     input: # "internal_output/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
#         windows = lambda wildcards: expand(
#             "internal_output/enriched_windows/{libname}.{{clip_sample_label}}.enriched_windows.tsv.gz", 
#             libname = manifest.loc[manifest['experiment']==wildcards.experiment, 'libname'].tolist()) # find all library of the same experiment
#     output:
#         reproducible_windows = "output/reproducible_enriched_windows/{experiment}.{clip_sample_label}.reproducible_enriched_windows.tsv.gz",
#         linear_bar = "output/figures/reproducible_enriched_windows/{experiment}.{clip_sample_label}.reproducible_enriched_window_counts.linear.pdf",
#         log_bar = "output/figures/reproducible_enriched_windows/{experiment}.{clip_sample_label}.reproducible_enriched_window_counts.log10.pdf"
#     params:
#         error_file = "stderr/{experiment}.{clip_sample_label}.find_reproducible_enriched_windows.err",
#         out_file = "stdout/{experiment}.{clip_sample_label}.find_reproducible_enriched_windows.out",
#         run_time = "5:00",
#         memory = "2000",
#     benchmark: "benchmarks/find_reproducible_enriched_windows/{experiment}.{clip_sample_label}.all_replicates.reproducible.txt"
#     container:
        # "docker://howardxu520/skipper:R_4.1.3_1"
        # shell:
#         "Rscript --vanilla {SCRIPT_PATH}/identify_reproducible_windows.R internal_output/enriched_windows/ {wildcards.clip_sample_label} " + (BLACKLIST if BLACKLIST is not None else "") 
