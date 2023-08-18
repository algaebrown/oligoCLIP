SCRIPT_PATH = config['SCRIPT_PATH']
FEATURE_ANNOTATIONS = config['FEATURE_ANNOTATIONS']
R_EXE = config['R_EXE']

rule sample_background_windows_by_region:
    input:
        enriched_windows = lambda wildcards: f"{wildcards.root_dir}/enriched_windows/{wildcards.libname}.{wildcards.sample_label}.enriched_windows.tsv.gz" if 'skipper' 
            in wildcards.root_dir else f"{wildcards.root_dir}/{wildcards.libname}.{wildcards.sample_label}.enriched_windows.tsv",
        all_windows = FEATURE_ANNOTATIONS,
    output:
        variable_windows = "{root_dir}/homer/region_matched_background/variable/{libname}.{sample_label}.sampled_variable_windows.bed.gz",
        fixed_windows = "{root_dir}/homer/region_matched_background/fixed/{libname}.{sample_label}.sampled_fixed_windows.bed.gz"
    params:
        error_out_file = "error_files/sample_background_windows_by_region.{libname}.{sample_label}.{root_dir}.err",
        out_file = "stdout/sample_background_windows_by_region.{libname}.{sample_label}.{root_dir}.out",
        run_time = "10:00",
        cores = 1,
        window_size =  75,
        outdir = lambda wildcards, output: output.variable_windows.split('variable')[0]
    benchmark: "benchmarks/sample_background_windows_by_region/{libname}.{sample_label}.{root_dir}.sample_background_windows_by_region.txt"
    shell:
        """
        {R_EXE} --vanilla {SCRIPT_PATH}/sample_matched_background_by_region.R \
            {input.enriched_windows} \
            {input.all_windows} \
            {params.window_size} \
            {params.outdir} \
            {wildcards.libname}.{wildcards.sample_label};
        """
rule run_homer:
    input:
        finemapped_windows = "{root_dir}/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
        background = "{root_dir}/homer/region_matched_background/fixed/{libname}.{sample_label}.sampled_fixed_windows.bed.gz",
        genome = config['GENOMEFA']
    output:
        report = "{root_dir}/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
        motif = "{root_dir}/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerMotifs.all.motifs"
    params:
        error_out_file = "error_files/run_homer.{libname}.{sample_label}.{signal_type}.{root_dir}.err",
        out_file = "stdout/run_homer.{libname}.{sample_label}.{signal_type}.{root_dir}.out",
        run_time = "40:00",
        cores = 1,
        outdir = lambda wildcards, output: output.motif.split('finemapped_results')[0]
    conda:
         "envs/homer.yaml"
    benchmark: "benchmarks/run_homer/{libname}.{sample_label}.{signal_type}.{root_dir}.txt"
    shell:
        """
        findMotifsGenome.pl <(less {input.finemapped_windows} | awk -v OFS=\"\t\" '{{print $4 \":\"$9,$1,$2+1,$3,$6}}') \
            {input.genome} {params.outdir}/finemapped_results/{wildcards.signal_type}/{wildcards.libname}.{wildcards.sample_label} \
            -preparsedDir {params.outdir}/preparsed -size given -rna -nofacts -S 20 -len 5,6,7,8,9 -nlen 1 \
            -bg <(zcat {input.background} | awk -v OFS=\"\t\" '{{print $4,$1,$2+1,$3,$6}}') 
        """

# use rule sample_background_windows_by_region as sample_background_windows_by_region_beta with:
#     input:
#         enriched_windows = "beta-mixture_CC/{libname}.{sample_label}.enriched_windows.tsv",
#         all_windows = FEATURE_ANNOTATIONS,
#     output:
#         variable_windows = "beta-mixture_CC/homer/region_matched_background/variable/{libname}.{sample_label}.sampled_variable_windows.bed.gz",
#         fixed_windows = "beta-mixture_CC/homer/region_matched_background/fixed/{libname}.{sample_label}.sampled_fixed_windows.bed.gz"
    

# use rule run_homer as run_homer_beta with:
#     input:
#         finemapped_windows = "beta-mixture_CC/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
#         background = "beta-mixture_CC/homer/region_matched_background/fixed/{libname}.{sample_label}.sampled_fixed_windows.bed.gz",
#         genome = config['GENOMEFA']
#     output:
#         report = "beta-mixture_CC/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
#         motif = "beta-mixture_CC/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerMotifs.all.motifs"
    

# use rule sample_background_windows_by_region as sample_background_windows_by_region_dmm with:
#     input:
#         enriched_windows = "DMM/{libname}.{sample_label}.enriched_windows.tsv",
#         all_windows = FEATURE_ANNOTATIONS,
#     output:
#         variable_windows = "DMM/homer/region_matched_background/variable/{libname}.{sample_label}.sampled_variable_windows.bed.gz",
#         fixed_windows = "DMM/homer/region_matched_background/fixed/{libname}.{sample_label}.sampled_fixed_windows.bed.gz"

# use rule run_homer as run_homer_dmm with:
#     input:
#         finemapped_windows = "DMM/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
#         background = "DMM/homer/region_matched_background/fixed/{libname}.{sample_label}.sampled_fixed_windows.bed.gz",
#         genome = config['GENOMEFA']
#     output:
#         report = "DMM/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
#         motif = "DMM/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerMotifs.all.motifs"
    

