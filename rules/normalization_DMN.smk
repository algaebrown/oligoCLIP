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

rule make_window_by_barcode_table:
    input:
        counts = expand("output/counts/genome/vectors/{libname}.{sample_label}.counts",
            libname = ["{libname}"],
            sample_label = rbps),
    output:
        counts = "DMN/table/{libname}.tsv.gz",
    params:
        error_out_file = "error_files/{libname}.window_by_barcode_table.err",
        out_file = "stdout/{libname}.window_by_barcode_table.out",
        run_time = "20:00",
        cores = 1
    shell:
        """
        paste -d '\t' {input.counts} | gzip  > {output.counts}
        """

rule make_read_count_summary:
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        counts = "DMN/table/{libname}.tsv.gz"
    output:
        region_summary =  "QC/read_count/{libname}.region.csv",
        type_summary =  "QC/read_count/{libname}.genetype.csv",
        name_summary =  "QC/read_count/{libname}.genename.csv",
        dist = "QC/read_count/{libname}.cosine_similarity.csv"
    params:
        error_out_file = "error_files/{libname}.read_count_summary.err",
        out_file = "stdout/{libname}.read_count_summary.out",
        run_time = "1:20:00",
        cores = 1
    run:
        import os
        print(output.region_summary)
        try:
            os.mkdir('QC/read_count')
        except Exception as e:
            print(e)
        import pandas as pd
        cnt = pd.read_csv(input.counts, sep = '\t')
        feature_annotations = pd.read_csv(input.feature_annotations, sep = '\t')

        df = pd.concat([feature_annotations, cnt], axis = 1)

        by_type = df.groupby(by = 'feature_type_top')[cnt.columns].sum()
        by_gene = df.groupby(by = 'gene_type_top')[cnt.columns].sum()
        by_name = df.groupby(by = 'gene_name')[cnt.columns].sum()

        by_type.to_csv(output.region_summary)
        by_gene.to_csv(output.type_summary)
        by_name.to_csv(output.name_summary)

        # distance
        from scipy.spatial.distance import pdist, squareform
        cov_filter = 10
        dist = squareform(pdist(cnt.loc[cnt.sum(axis = 1)>cov_filter].T, 'cosine'))

        dist_df = pd.DataFrame(1-dist, columns = cnt.columns, index = cnt.columns)
        dist_df.to_csv(output.dist)
