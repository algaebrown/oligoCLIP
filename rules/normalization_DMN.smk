import pandas as pd
R_EXE = config['R_EXE']
BLACKLIST = config['BLACKLIST'] if 'BLACKLIST' in config else None
rbps = config['rbps']
experiments = config['experiments']
libnames = config['libnames']
SCRIPT_PATH=config['SCRIPT_PATH']
manifest = pd.read_table(config['MANIFEST'], index_col = False, sep = ',')
UNINFORMATIVE_READ = 3 - int(config['INFORMATIVE_READ']) # whether read 1 or read 2 is informative



def libname_to_experiment(libname):
    return manifest.loc[manifest['libname']==libname, 'experiment'].iloc[0]
def experiment_to_libname(experiment):
    libnames = manifest.loc[manifest['experiment']==experiment, 'libname'].tolist()
    assert len(libnames)>0
    return libnames


rule fit_beta_mixture_model_internal:
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        table = lambda wildcards: f"internal_output/counts/genome/bgtables/internal/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
    output:
        "internal_output/DMN/{libname}.{clip_sample_label}.fit.rda",
        "internal_output/DMN/{libname}.{clip_sample_label}.goodness_of_fit.pdf",
        "internal_output/DMN/{libname}.{clip_sample_label}.weights.tsv",
        "internal_output/DMN/{libname}.{clip_sample_label}.alpha.tsv",
        "internal_output/DMN/{libname}.{clip_sample_label}.null.alpha.tsv",
        "internal_output/DMN/{libname}.{clip_sample_label}.mixture_weight.tsv",
    params:
        error_out_file = "error_files/{libname}.{clip_sample_label}.internal.DMN.err",
        out_file = "stdout/{libname}.{clip_sample_label}.internal.DMN.err",
        run_time = "00:40:00",
        memory = "10000",
        job_name = "DMN_internal",
        cores = "1",
    conda:
        "envs/DMN.yaml"
    benchmark: "benchmarks/DMN/fit_model.{libname}.{clip_sample_label}"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/fit_DMN.R \
            {input.table} \
            {input.feature_annotations} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {wildcards.libname}.internal \
            internal_output/DMN/ \
            {wildcards.libname}.{wildcards.clip_sample_label} &> {params.out_file}
        """

rule analyze_beta_mixture_results:
    input:
        "internal_output/DMN/{libname}.{clip_sample_label}.fit.rda",
        "internal_output/DMN/{libname}.{clip_sample_label}.goodness_of_fit.pdf",
        "internal_output/DMN/{libname}.{clip_sample_label}.weights.tsv",
        "internal_output/DMN/{libname}.{clip_sample_label}.alpha.tsv",
        "internal_output/DMN/{libname}.{clip_sample_label}.null.alpha.tsv",
        "internal_output/DMN/{libname}.{clip_sample_label}.mixture_weight.tsv",
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        table = lambda wildcards: f"internal_output/counts/genome/bgtables/internal/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
    output:
        "internal_output/DMN/{libname}.{clip_sample_label}.enriched_window.tsv"
    params:
        error_out_file = "error_files/{libname}.{clip_sample_label}.analyze.DMN.err",
        out_file = "stdout/{libname}.{clip_sample_label}.analyze.DMN.err",
        run_time = "00:40:00",
        memory = "10000",
        job_name = "DMN_analysis",
        cores = "1",
    conda:
        "envs/metadensity.yaml"
    benchmark: "benchmarks/DMN/analyze.{libname}.{clip_sample_label}"
    shell:
        """
        python {SCRIPT_PATH}/analyze_betabinom_mixture.py \
            internal_output/DMN/ \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {input.table} \
            {input.feature_annotations}
        """

rule analyze_beta_mixture_results_most_enrich:
    input:
        "internal_output/DMN/{libname}.{clip_sample_label}.fit.rda",
        "internal_output/DMN/{libname}.{clip_sample_label}.goodness_of_fit.pdf",
        "internal_output/DMN/{libname}.{clip_sample_label}.weights.tsv",
        "internal_output/DMN/{libname}.{clip_sample_label}.alpha.tsv",
        "internal_output/DMN/{libname}.{clip_sample_label}.null.alpha.tsv",
        "internal_output/DMN/{libname}.{clip_sample_label}.mixture_weight.tsv",
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        table = lambda wildcards: f"internal_output/counts/genome/bgtables/internal/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
    output:
        "internal_output/DMN/most_enriched_selected/{libname}.{clip_sample_label}.enriched_window.tsv"
    params:
        error_out_file = "error_files/{libname}.{clip_sample_label}.analyze.DMN.err",
        out_file = "stdout/{libname}.{clip_sample_label}.analyze.DMN.err",
        run_time = "00:40:00",
        memory = "10000",
        job_name = "DMN_analysis",
        cores = "1",
    conda:
        "envs/metadensity.yaml"
    benchmark: "benchmarks/DMN/analyze.{libname}.{clip_sample_label}"
    shell:
        """
        python {SCRIPT_PATH}/analyze_betabinom_mixture_most_enriched.py \
            internal_output/DMN/ \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {input.table} \
            {input.feature_annotations}
        """

rule fit_beta_mixture_model_another_lib:
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        table = lambda wildcards: f"output/counts/genome/bgtables/{wildcards.bg_sample_label}/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz"
    output:
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.fit.rda",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.goodness_of_fit.pdf",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.weights.tsv",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.alpha.tsv",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.null.alpha.tsv",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.mixture_weight.tsv",
    params:
        error_out_file = "error_files/{libname}.{clip_sample_label}.{bg_sample_label}.internal.DMN.err",
        out_file = "stdout/{libname}.{clip_sample_label}.{bg_sample_label}.internal.DMN.err",
        run_time = "00:40:00",
        memory = "10000",
        job_name = "DMN_internal",
        cores = "1",
    conda:
        "envs/DMN.yaml"
    benchmark: "benchmarks/DMN/fit_model.{libname}.{clip_sample_label}.{bg_sample_label}"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/fit_DMN.R \
            {input.table} \
            {input.feature_annotations} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {wildcards.libname}.{wildcards.bg_sample_label} \
            output/DMN/{wildcards.bg_sample_label} \
            {wildcards.libname}.{wildcards.clip_sample_label} 
        """


rule analyze_beta_mixture_results_another_lib:
    input:
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.fit.rda",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.goodness_of_fit.pdf",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.weights.tsv",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.alpha.tsv",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.null.alpha.tsv",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.mixture_weight.tsv",
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        table = lambda wildcards: f"output/counts/genome/bgtables/{wildcards.bg_sample_label}/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz"
    output:
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.enriched_window.tsv"
    params:
        error_out_file = "error_files/{libname}.{clip_sample_label}.{bg_sample_label}.analyze.DMN.err",
        out_file = "stdout/{libname}.{clip_sample_label}.{bg_sample_label}.analyze.DMN.err",
        run_time = "00:40:00",
        memory = "10000",
        job_name = "DMN_analysis",
        cores = "1",
    conda:
        "envs/metadensity.yaml"
    benchmark: "benchmarks/DMN/analyze.{libname}.{clip_sample_label}.{bg_sample_label}"
    shell:
        """
        python {SCRIPT_PATH}/analyze_betabinom_mixture.py \
            output/DMN/{wildcards.bg_sample_label} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {input.table} \
            {input.feature_annotations}
        """

rule analyze_beta_mixture_results_another_lib_most_enrich:
    input:
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.fit.rda",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.goodness_of_fit.pdf",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.weights.tsv",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.alpha.tsv",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.null.alpha.tsv",
        "output/DMN/{bg_sample_label}/{libname}.{clip_sample_label}.mixture_weight.tsv",
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        table = lambda wildcards: f"output/counts/genome/bgtables/{wildcards.bg_sample_label}/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz"
    output:
        "output/DMN/{bg_sample_label}/most_enriched_selected/{libname}.{clip_sample_label}.enriched_window.tsv"
    params:
        error_out_file = "error_files/{libname}.{clip_sample_label}.{bg_sample_label}.analyze.DMN.err",
        out_file = "stdout/{libname}.{clip_sample_label}.{bg_sample_label}.analyze.DMN.err",
        run_time = "00:40:00",
        memory = "10000",
        job_name = "DMN_analysis",
        cores = "1",
    conda:
        "envs/metadensity.yaml"
    benchmark: "benchmarks/DMN/analyze.{libname}.{clip_sample_label}.{bg_sample_label}"
    shell:
        """
        python {SCRIPT_PATH}/analyze_betabinom_mixture_most_enriched.py \
            output/DMN/{wildcards.bg_sample_label} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {input.table} \
            {input.feature_annotations}
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
