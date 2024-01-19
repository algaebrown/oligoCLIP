import pandas as pd
locals().update(config)

manifest = pd.read_table(config['MANIFEST'], index_col = False, sep = ',')



def libname_to_experiment(libname):
    return manifest.loc[manifest['libname']==libname, 'experiment'].iloc[0]
def experiment_to_libname(experiment):
    libnames = manifest.loc[manifest['experiment']==experiment, 'libname'].tolist()
    assert len(libnames)>0
    return libnames


rule fit_beta_mixture_model_CC:
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        table = lambda wildcards: "counts_CC/genome/bgtables/internal/"+libname_to_experiment(wildcards.libname)+'.'+wildcards.clip_sample_label+".tsv.gz",
    output:
        "beta-mixture_CC/{libname}.{clip_sample_label}.fit.rda",
        "beta-mixture_CC/{libname}.{clip_sample_label}.goodness_of_fit.pdf",
        "beta-mixture_CC/{libname}.{clip_sample_label}.weights.tsv",
        "beta-mixture_CC/{libname}.{clip_sample_label}.alpha.tsv",
        "beta-mixture_CC/{libname}.{clip_sample_label}.null.alpha.tsv",
        "beta-mixture_CC/{libname}.{clip_sample_label}.mixture_weight.tsv",
    params:
        error_out_file = "error_files/fit_betaCC.{libname}.{clip_sample_label}.err",
        out_file = "stdout/fit_betaCC.{libname}.{clip_sample_label}.out",
        run_time = "02:40:00",
        cores = "1",
        memory = 80000,
        root_folder = lambda wildcards, output: Path(output[0]).parent
    conda:
        "envs/DMN.yaml"
    benchmark: "benchmarks/DMM/fit_beta.{libname}.{clip_sample_label}"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/fit_DMN.R \
            {input.table} \
            {input.feature_annotations} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {wildcards.libname}.internal \
            {params.root_folder} \
            {wildcards.libname}.{wildcards.clip_sample_label} 
        """

rule analyze_beta_mixture_results_CC:
    input:
        "beta-mixture_CC/{libname}.{clip_sample_label}.fit.rda",
        "beta-mixture_CC/{libname}.{clip_sample_label}.goodness_of_fit.pdf",
        "beta-mixture_CC/{libname}.{clip_sample_label}.weights.tsv",
        "beta-mixture_CC/{libname}.{clip_sample_label}.alpha.tsv",
        "beta-mixture_CC/{libname}.{clip_sample_label}.null.alpha.tsv",
        "beta-mixture_CC/{libname}.{clip_sample_label}.mixture_weight.tsv",
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        table = lambda wildcards: "counts_CC/genome/bgtables/internal/"+libname_to_experiment(wildcards.libname)+'.'+wildcards.clip_sample_label+".tsv.gz",
    output:
        "beta-mixture_CC/{libname}.{clip_sample_label}.enriched_windows.tsv"
    params:
        error_out_file = "error_files/analyze_beta_CC.{libname}.{clip_sample_label}.err",
        out_file = "stdout/analyze_beta_CC.{libname}.{clip_sample_label}.err",
        run_time = "00:40:00",
        cores = "1",
        root_folder = lambda wildcards, output: Path(output[0]).parent,
        memory = 40000,
    conda:
        "envs/tensorflow.yaml"
    benchmark: "benchmarks/DMN/analyze.{libname}.{clip_sample_label}"
    shell:
        """
        python {SCRIPT_PATH}/analyze_betabinom_mixture_most_enriched.py \
            {params.root_folder} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {input.table} \
            {input.feature_annotations}
        """

rule fit_beta_mixture_model_another_lib:
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        table = lambda wildcards: f"counts/genome/bgtables/{wildcards.bg_sample_label}/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz"
    output:
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.fit.rda",
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.goodness_of_fit.pdf",
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.weights.tsv",
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.alpha.tsv",
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.null.alpha.tsv",
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.mixture_weight.tsv",
    params:
        error_out_file = "error_files/fit_beta.{bg_sample_label}.{libname}.{clip_sample_label}.err",
        out_file = "error_files/fit_beta.{bg_sample_label}.{libname}.{clip_sample_label}.out",
        run_time = "00:40:00",
        cores = "1",
        root_folder = lambda wildcards, output: Path(output[0]).parent,
        memory = 80000,
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
            {params.root_folder} \
            {wildcards.libname}.{wildcards.clip_sample_label} 
        """


rule analyze_beta_mixture_results_another_lib:
    input:
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.fit.rda",
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.goodness_of_fit.pdf",
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.weights.tsv",
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.alpha.tsv",
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.null.alpha.tsv",
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.mixture_weight.tsv",
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        table = lambda wildcards: f"counts/genome/bgtables/{wildcards.bg_sample_label}/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz"
    output:
        "beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.enriched_windows.tsv"
    params:
        error_out_file = "error_files/analyze_beta.{bg_sample_label}.{libname}.{clip_sample_label}.err",
        out_file = "stdout/analyze_beta.{bg_sample_label}.{libname}.{clip_sample_label}.err",
        run_time = "00:40:00",
        cores = "1",
        root_folder = lambda wildcards, output: Path(output[0]).parent,
        memory = 40000,
    conda:
        "envs/tensorflow.yaml"
    benchmark: "benchmarks/DMN/analyze.{libname}.{clip_sample_label}.{bg_sample_label}"
    shell:
        """
        python {SCRIPT_PATH}/analyze_betabinom_mixture_most_enriched.py \
            {params.root_folder} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {input.table} \
            {input.feature_annotations}
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
#     shell:
#         "Rscript --vanilla {SCRIPT_PATH}/identify_reproducible_windows.R internal_output/enriched_windows/ {wildcards.clip_sample_label} " + (BLACKLIST if BLACKLIST is not None else "") 

rule make_window_by_barcode_table:
    input:
        counts = expand("counts/genome/vectors/{libname}.{sample_label}.counts",
            libname = ["{libname}"],
            sample_label = rbps),
    output:
        counts = "counts/genome/megatables/{libname}.tsv.gz",
    params:
        error_out_file = "error_files/window_by_barcode_table.{libname}.err",
        out_file = "stdout/window_by_barcode_table.{libname}.out",
        run_time = "20:00",
        cores = 1,
        memory = 40000,
    shell:
        """
        paste -d '\t' {input.counts} | gzip  > {output.counts}
        """

rule fit_DMN:
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        table = "counts/genome/megatables/{libname}.tsv.gz",
    output:
        "DMM/{libname}.goodness_of_fit.pdf",
        "DMM/{libname}.alpha.tsv",
        "DMM/{libname}.null.alpha.tsv",
        "DMM/{libname}.mixture_weight.tsv",
        "DMM/{libname}.weights.tsv"
    params:
        error_out_file = "error_files/fit_DMM.{libname}.err",
        out_file = "stdout/DMM.{libname}.internal.out",
        run_time = "48:00:00",
        memory = 40000,
        cores = "4",
        root_folder = "DMM",
    benchmark: "benchmarks/DMM/fit.{libname}"
    conda:
        "envs/DMN.yaml"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/fit_DMN_multidimen.R \
            {input.table} \
            {input.feature_annotations} \
            {params.root_folder} \
            {wildcards.libname} 
        """

rule analyze_DMN:
    input:
        "DMM/{libname}.alpha.tsv",
        "DMM/{libname}.mixture_weight.tsv",
        "counts/genome/megatables/{libname}.tsv.gz",
        'QC/mapping_stats.csv',
        "DMM/{libname}.weights.tsv",
        'mask/{libname}.genome_mask.csv'
    output:
        expand("DMM/{libname}.{sample_labels}.enriched_windows.tsv", sample_labels = rbps, libname = ["{libname}"]),
        "DMM/{libname}.megaoutputs.tsv"
    params:
        error_out_file = "error_files/analyze_DMN.{libname}.err",
        out_file = "stdout/analyze_DMN.{libname}.out",
        run_time = "2:00:00",
        cores = "1",
        memory = 80000,
    benchmark: "benchmarks/DMM/analysis.{libname}"
    conda:
        "envs/tensorflow.yaml"
    shell:
        """
        python {SCRIPT_PATH}/analyze_DMM.py {wildcards.libname}
        """

rule softmask:
    input:
        "counts/genome/megatables/{libname}.tsv.gz",
        "counts/repeats/megatables/name/{libname}.tsv.gz",
        rep_annotation = config['REPEAT_TABLE'],
        genomic_annotation = config['FEATURE_ANNOTATIONS'],
    output:
        genomic_dev_zscore = 'mask/{libname}.genome_deviation_zscore.csv',
        genome_mask = 'mask/{libname}.genome_mask.csv',
        repeat_dev_zscore = 'mask/{libname}.repeat_deviation_zscore.csv',
        repeat_mask = 'mask/{libname}.repeat_mask.csv', #True means zscore > 2
    params:
        error_out_file = "error_files/softmask.{libname}.err",
        out_file = "stdout/softmask.{libname}.out",
        run_time = "2:00:00",
        cores = "1",
        memory = 80000,
    benchmark: "benchmarks/DMM/softmask.{libname}"
    conda:
        "envs/tensorflow.yaml"
    shell:
        """
        python {SCRIPT_PATH}/softmask_noisy_region.py . {wildcards.libname} {input.genomic_annotation} {input.rep_annotation}
        """


###### GC-aware rules for external normalization, calling peaks when comparing libraries amplifed seperately #####
rule fit_beta_gc_aware:
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        table = lambda wildcards: f'counts_external/genome/{wildcards.external_label}/'+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
        gc = config['PARTITION'].replace('bed.gz', 'nuc.gz')
    output:
        expand("beta-mixture_external/{external_label}/{libname}.{clip_sample_label}.gc{index}.goodness_of_fit.pdf",
            external_label = ['{external_label}'], libname = ['{libname}'], clip_sample_label = ['{clip_sample_label}'],
            index = list(range(1,11))),
        expand("beta-mixture_external/{external_label}/{libname}.{clip_sample_label}.gc{index}.alpha.tsv",
            external_label = ['{external_label}'], libname = ['{libname}'], clip_sample_label = ['{clip_sample_label}'],
                index = list(range(1,11))),
        expand("beta-mixture_external/{external_label}/{libname}.{clip_sample_label}.gc{index}.null.alpha.tsv",
            external_label = ['{external_label}'], libname = ['{libname}'], clip_sample_label = ['{clip_sample_label}'],
                index = list(range(1,11))),
        expand("beta-mixture_external/{external_label}/{libname}.{clip_sample_label}.gc{index}.mixture_weight.tsv",
            external_label = ['{external_label}'], libname = ['{libname}'], clip_sample_label = ['{clip_sample_label}'],
                index = list(range(1,11))),
        expand("beta-mixture_external/{external_label}/{libname}.{clip_sample_label}.gc{index}.weights.tsv",
            external_label = ['{external_label}'], libname = ['{libname}'], clip_sample_label = ['{clip_sample_label}'],
                index = list(range(1,11)))
    params:
        error_out_file = "error_files/fit_beta_gcaware.{libname}.{clip_sample_label}.{external_label}.err",
        out_file = "stdout/beta-mixture_GC.{libname}.{clip_sample_label}.{external_label}.internal.out",
        run_time = "2:00:00",
        cores = "2",
        memory = 80000,
        root_folder = "beta-mixture_external/{external_label}",

    benchmark: "benchmarks/DMM/fit-beta-mixture_GC.{libname}.{clip_sample_label}.{external_label}"
    conda:
        "envs/DMN.yaml"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/fit_DMN_gcaware.R \
            {input.table} \
            {input.feature_annotations} \
            {input.gc} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            external.{wildcards.external_label} \
            {params.root_folder} \
            {wildcards.libname}.{wildcards.clip_sample_label}
        """

rule analyze_beta_GC_aware:
    input:
        rules.fit_beta_gc_aware.output,
        table = lambda wildcards: f"counts_external/genome/{wildcards.external_label}/"+libname_to_experiment(wildcards.libname)+f".{wildcards.clip_sample_label}.tsv.gz",
        feature_annotations = config['FEATURE_ANNOTATIONS'],
    output:
        "beta-mixture_external/{external_label}/{libname}.{clip_sample_label}.enriched_windows.tsv",
        "beta-mixture_external/{external_label}/{libname}.{clip_sample_label}.window_score.tsv",
    params:
        error_out_file = "error_files/analyze_beta_gcaware.{libname}.{clip_sample_label}.{external_label}.err",
        out_file = "stdout/analyze_beta-mixture_GC.{libname}.{clip_sample_label}.{external_label}.out",
        run_time = "1:00:00",
        cores = "1",
        root_folder = "beta-mixture_external/{external_label}",
        memory = 80000,
    benchmark: "benchmarks/DMM/analyze-beta-mixture_GC.{libname}.{clip_sample_label}.{external_label}"
    conda:
        "envs/tensorflow.yaml"
    shell:
        """
        python {SCRIPT_PATH}/analyze_betabinom_mixture_most_enriched_gcaware.py \
            {params.root_folder} \
            {wildcards.libname}.{wildcards.clip_sample_label} \
            {input.table} \
            {input.feature_annotations} \
            external.{wildcards.external_label}
        """