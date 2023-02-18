# zcat ../repeatmasker.grch38.sor6t.unique.bed.gz | paste - output/counts/repeats/vectors/
CHROM_SIZES = config['CHROM_SIZES']
PARTITION = config['PARTITION']
GENOME = config['GENOMEFA']
FEATURE_ANNOTATIONS = config['FEATURE_ANNOTATIONS']
ACCESSION_RANKINGS = config['ACCESSION_RANKINGS']
UNINFORMATIVE_READ = 3 - int(config['INFORMATIVE_READ']) # whether read 1 or read 2 is informative

rule partition_bam_reads:
    input:
        chrom_size = config['CHROM_SIZES'],
        bam = "output/bams/dedup/genome/{replicate_label}.genome.Aligned.sort.dedup.bam",        
        region_partition = PARTITION,
    output:
        counts = "output/counts/genome/vectors/{replicate_label}.counts",
    params:
        error_out_file = "stderr/{replicate_label}.partition_bam_reads.err",
        out_file = "stdout/{replicate_label}.partition_bam_reads.out",
        run_time = "20:00",
        cores = "1",
        memory = "10000",
        job_name = "partition_bam_reads",
        replicate_label = "{replicate_label}"
    benchmark: "benchmarks/counts/unassigned_experiment.{replicate_label}.partition_bam_reads.txt"
    shell:
        "bedtools bamtobed -i {input.bam} | awk '($1 != \"chrEBV\") && ($4 !~ \"/{UNINFORMATIVE_READ}$\")' | bedtools flank -s -l 1 -r 0 -g {input.chrom_size} -i - | bedtools shift -p -1 -m 1 -g {input.chrom_size} -i - | bedtools coverage -counts -s -a {input.region_partition} -b - | cut -f 7 | awk 'BEGIN {{print \"{params.replicate_label}\"}} {{print}}' > {output.counts};"

rule make_genome_count_table:
    input:
        partition=config['PARTITION'],
        replicate_counts = lambda wildcards: expand("output/counts/genome/vectors/{replicate_label}.counts", replicate_label = experiment_to_replicate_labels[wildcards.experiment_label]),
    output:
        count_table = "output/counts/genome/tables/{experiment_label}.tsv.gz",
    params:
        error_out_file = "stderr/{experiment_label}.make_count_table.err",
        out_file = "stdout/{experiment_label}.make_count_table.out",
        run_time = "00:05:00",
        cores = "1",
        memory = "200",
        job_name = "make_genome_count_table"
    benchmark: "benchmarks/counts/{experiment_label}.all_replicates.make_genome_count_table.txt"
    shell:
        "paste <(zcat {input.partition} | awk -v OFS=\"\\t\" 'BEGIN {{print \"chr\\tstart\\tend\\tname\\tscore\\tstrand\"}} {{print $1,$2,$3,$4,$5,$6}}' ) {input.replicate_counts} | gzip -c > {output.count_table}"
         #"paste                          <(awk -v OFS=\"\\t\" 'BEGIN {{print \"chr\\tstart\\tend\\tname\\tscore\\tstrand\"}} {{print}}' {input.partition})                 {input.replicate_counts} | gzip -c > {output.count_table};"

rule calc_partition_nuc:
    input:
        partition = PARTITION,
        genome = GENOME
    output:
        nuc = PARTITION.replace(".bed", ".nuc")
    params:
        error_out_file = "stderr/calc_partition_nuc.err",
        out_file = "stdout/calc_partition_nuc.out",
        run_time = "00:10:00",
        memory = "1000",
        job_name = "calc_partition_nuc"
    benchmark: "benchmarks/partition_nuc.txt"
    shell:
        "bedtools nuc -s -fi {input.genome} -bed {input.partition} > {output.nuc}"

rule fit_input_betabinomial_model:
    input:
        nuc = PARTITION.replace(".bed", ".nuc"),
        table = "output/counts/genome/tables/{experiment_label}.tsv.gz"
    output:
        coef = "output/input_model_coef/{experiment_label}.{input_replicate_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/input_distributions/{{experiment_label}}.{{input_replicate_label}}.{other_label}.input_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment_label][wildcards.Input_replicate_label])
    params:
        error_out_file = "stderr/{experiment_label}_{input_replicate_label}.fit_input_betabinom.err",
        out_file = "stdout/{experiment_label}_{input_replicate_label}.fit_input_betabinom.out",
        run_time = "1:00:00",
        memory = "1000",
        job_name = "fit_input_betabinomial_model"
    benchmark: "benchmarks/betabinomial/{experiment_label}.{input_replicate_label}.fit_input.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/fit_input_betabinom.R {input.nuc} {input.table} {wildcards.experiment_label} {wildcards.input_replicate_label}"

rule fit_clip_betabinomial_model:
    input:
        nuc = PARTITION.replace(".bed", ".nuc"),
        table = "output/counts/genome/tables/{experiment_label}.tsv.gz"
    output:
        coef = "output/clip_model_coef/{experiment_label}.{clip_replicate_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/clip_distributions/{{experiment_label}}.{{clip_replicate_label}}.{other_label}.clip_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment_label][wildcards.Input_replicate_label])
    params:
        error_out_file = "stderr/{experiment_label}_{clip_replicate_label}.fit_clip_betabinomial_model.err",
        out_file = "stdout/{experiment_label}_{clip_replicate_label}.fit_clip_betabinomial_model.out",
        run_time = "1:00:00",
        memory = "1000",
        job_name = "fit_clip_betabinomial_model"
    benchmark: "benchmarks/fit_clip_betabinomial_model/{experiment_label}.{clip_replicate_label}.fit_clip.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/fit_clip_betabinom.R {input.nuc} {input.table} {wildcards.experiment_label} {wildcards.clip_replicate_label}"

rule call_enriched_windows:
    input:
        feature_annotations = FEATURE_ANNOTATIONS,
        accession_rankings = ACCESSION_RANKINGS,
        nuc = PARTITION.replace(".bed", ".nuc"),
        table = "output/counts/genome/tables/{experiment_label}.tsv.gz",
        parameters = lambda wildcards: "output/input_model_coef/{experiment_label}." + clip_to_input_replicate_label[wildcards.clip_replicate_label] + ".tsv",
    output:
        "output/threshold_scan/{experiment_label}.{clip_replicate_label}.threshold_data.tsv",
        "output/tested_windows/{experiment_label}.{clip_replicate_label}.tested_windows.tsv.gz",
        "output/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_windows.tsv.gz",
        "output/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_feature_data.tsv",
        "output/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_transcript_data.tsv",
        "output/enrichment_summaries/{experiment_label}.{clip_replicate_label}.enriched_window_gene_data.tsv",
        "output/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_fractions_feature_data.tsv",
        "output/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds_feature_data.tsv",
        "output/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds_transcript_data.tsv",
        "output/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds_feature_gc_data.tsv",
        "output/figures/threshold_scan/{experiment_label}.{clip_replicate_label}.threshold_scan.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_coverage.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_rates.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_counts.linear.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_counts.log10.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_odds.feature.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_odds.all_transcript_types.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_odds.select_transcript_types.pdf",
        "output/figures/enriched_windows/{experiment_label}.{clip_replicate_label}.enriched_window_counts.per_gene_feature.pdf",
        "output/figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_fractions.feature.pdf",
        "output/figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds.feature.pdf",
        "output/figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds.all_transcript_types.pdf",
        "output/figures/all_reads/{experiment_label}.{clip_replicate_label}.all_reads_odds.feature_gc.pdf"
    params:
        input_replicate_label = lambda wildcards: clip_to_input_replicate_label[wildcards.clip_replicate_label],
        error_out_file = "stderr/{experiment_label}_{clip_replicate_label}.call_enriched_windows.err",
        out_file = "stdout/{experiment_label}_{clip_replicate_label}.call_enriched_windows.out",
        run_time = "00:25:00",
        memory = "1000",
        job_name = "call_enriched_windows"
    benchmark: "benchmarks/call_enriched_windows/{experiment_label}.{clip_replicate_label}.call_enriched_windows.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/call_enriched_windows.R {input.nuc} {input.table} {input.accession_rankings} {input.feature_annotations} {input.parameters} {params.input_replicate_label} {wildcards.clip_replicate_label} {wildcards.experiment_label}.{wildcards.clip_replicate_label}"

rule check_window_concordance:
    input:
        windows = lambda wildcards: expand("output/tested_windows/{{experiment_label}}.{clip_replicate_label}.tested_windows.tsv.gz", clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label])
    output:
        "output/figures/enrichment_reproducibility/{experiment_label}.enrichment_reproducibility.pdf",
        "output/enrichment_reproducibility/{experiment_label}.enrichment_reproducibility.tsv"
    params:
        error_out_file = "stderr/{experiment_label}.check_window_concordance.err",
        out_file = "stdout/{experiment_label}.check_window_concordance.out",
        run_time = "0:15:00",
        memory = "1000",
        job_name = "check_window_concordance"
    benchmark: "benchmarks/check_window_concordance/{experiment_label}.all_replicates.concordance.txt"
    shell:
        "{R_EXE} --vanilla {TOOL_DIR}/check_window_concordance.R output/tested_windows {wildcards.experiment_label}"

rule find_reproducible_enriched_windows:
    input:
        windows = lambda wildcards: expand("output/enriched_windows/{{experiment_label}}.{clip_replicate_label}.enriched_windows.tsv.gz", clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label])
    output:
        reproducible_windows = "output/reproducible_enriched_windows/{experiment_label}.reproducible_enriched_windows.tsv.gz"
    params:
        error_out_file = "stderr/{experiment_label}.find_reproducible_enriched_windows.err",
        out_file = "stdout/{experiment_label}.find_reproducible_enriched_windows.out",
        run_time = "1:00:00",
        memory = "1000",
        job_name = "find_reproducible_enriched_windows"
    benchmark: "benchmarks/find_reproducible_enriched_windows/{experiment_label}.all_replicates.reproducible.txt"
    shell:
        "zcat {input.windows} | awk -v OFS=\"\\t\" "
        "'($4 in common_fields) {{reproduced[$4] = 1}} "
        "{{  input_sum[$4] = input_sum[$4] + $7; "
            "clip_sum[$4] = clip_sum[$4] + $8; "
            "or_max[$4] = or_max[$4] > $12 ? or_max[$4] : $12; "
            "p_max[$4] = p_max[$4] > $13 ? p_max[$4] : $13; "
            "q_max[$4] = q_max[$4] > $14 ? q_max[$4] : $14; "
            "common_fields[$4] = $1 \"\\t\" $2 \"\\t\" $3 \"\\t\" $4 \"\\t\" $5 \"\\t\" $6 \"\\t\" $9 \"\\t\" $10 \"\\t\" $11 "
                "\"\\t\" $15 \"\\t\" $16 \"\\t\" $17 \"\\t\" $18 \"\\t\" $19 \"\\t\" $20 \"\\t\" $21 \"\\t\" $22 \"\\t\" $23 \"\\t\" $24 \"\\t\" $25 }} "
        "END {{for(id in reproduced) {{print common_fields[id], input_sum[id], clip_sum[id], or_max[id], p_max[id], q_max[id]}} }}' | "
            "sort -nrk 12,12 | "
            "awk 'BEGIN {{print \"chr\\tstart\\tend\\tname\\tscore\\tstrand\\tpct_gc\\tgc_bin\\tbaseline_LOR\\t"
            "feature_id\\tfeature_bin\\tfeature_type_top\\tfeature_types\\tgene_name\\tgene_id\\ttranscript_ids\\tgene_type_top\\ttranscript_type_top\\tgene_types\\t"
            "input_sum\\tclip_sum\\tenrichment_LOR_max\\tp_max\\tq_max\"}} $1 != \"chr\" {{print}}' | "
        "gzip > {output.reproducible_windows}"


