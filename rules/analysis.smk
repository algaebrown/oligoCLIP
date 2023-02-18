SCRIPT_PATH = config['SCRIPT_PATH']
FEATURE_ANNOTATIONS = config['FEATURE_ANNOTATIONS']
R_EXE = config['R_EXE']
# rule make_homer_foreground_background:
#     input:
#         enriched_windows = "internal_output/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
#         background = "internal_output/tested_windows/{libname}.{clip_sample_label}.tested_windows.tsv.gz",
#         annotation = config['FEATURE_ANNOTATIONS']
#     output:
#         expand('homer/inputs/{libname}.{clip_sample_label}.{region}.{type}.homer',
#         libname = ['{libname}'],
#         clip_sample_label = '{clip_sample_label}',
#         region = ['CDS', 'UTR5', 'INTRON', 'UTR3', 'ALL'],
#         type = ['enrich', 'tested'])
#     params:
#         error_out_file = "stderr/{libname}.{clip_sample_label}.make_homer.err",
#         out_file = "stdout/{libname}.{clip_sample_label}.make_homer.out",
#         run_time = "00:40:00",
#         d_log_odds_cutoff = 2,
#         cores = 1
#     conda:
#         "envs/metadensity.yaml"
#     shell:
#         """
#         python {SCRIPT_PATH}/prepare_homer.py {input.enriched_windows} {input.background} {params.d_log_odds_cutoff} {input.annotation} homer/inputs/{wildcards.libname}.{wildcards.clip_sample_label}
#         """

rule sample_background_windows_by_region:
    input:
        enriched_windows = "internal_output/enriched_windows/{libname}.{sample_label}.enriched_windows.tsv.gz",
        all_windows = FEATURE_ANNOTATIONS,
    output:
        variable_windows = "internal_output/homer/region_matched_background/variable/{libname}.{sample_label}.sampled_variable_windows.bed.gz",
        fixed_windows = "internal_output/homer/region_matched_background/fixed/{libname}.{sample_label}.sampled_fixed_windows.bed.gz"
    params:
        error_out_file = "stderr/{libname}.{sample_label}.sample_background_windows_by_region.err",
        out_file = "stdout/{libname}.{sample_label}.sample_background_windows_by_region.out",
        run_time = "10:00",
        memory = "3000",
        cores = 1,
        job_name = "sample_background_windows"
    benchmark: "benchmarks/sample_background_windows_by_region/{libname}.{sample_label}.sample_background_windows_by_region.txt"
    shell:
        """
        {R_EXE} --vanilla {SCRIPT_PATH}/sample_matched_background_by_region.R \
            {input.enriched_windows} \
            {input.all_windows} \
            75 \
            internal_output/homer/region_matched_background \
            {wildcards.libname}.{wildcards.sample_label};
        """
rule run_homer:
    input:
        finemapped_windows = "internal_output/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
        background = "internal_output/homer/region_matched_background/fixed/{libname}.{sample_label}.sampled_fixed_windows.bed.gz",
        genome = config['GENOMEFA']
    output:
        report = "internal_output/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
        motif = "internal_output/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerMotifs.all.motifs"
    params:
        error_out_file = "stderr/{libname}.{sample_label}.{signal_type}.run_homer.err",
        out_file = "stdout/{libname}.{sample_label}.{signal_type}.run_homer.out",
        run_time = "40:00",
        memory = "2000",
        cores = 1,
        job_name = "run_homer"
    conda:
         "envs/homer.yaml"
    benchmark: "benchmarks/run_homer/{libname}.{sample_label}.{signal_type}.txt"
    shell:
        """
        findMotifsGenome.pl <(less {input.finemapped_windows} | awk -v OFS=\"\t\" '{{print $4 \":\"$9,$1,$2+1,$3,$6}}') \
            {input.genome} internal_output/homer/finemapped_results/{wildcards.signal_type}/{wildcards.libname}.{wildcards.sample_label} \
            -preparsedDir internal_output/homer/preparsed -size given -rna -nofacts -S 20 -len 5,6,7,8,9 -nlen 1 \
            -bg <(zcat {input.background} | awk -v OFS=\"\t\" '{{print $4,$1,$2+1,$3,$6}}') 
        """

rule sample_background_windows_by_region_dmn:
    input:
        enriched_windows = "internal_output/DMN/{libname}.{sample_label}.enriched_window.tsv",
        all_windows = FEATURE_ANNOTATIONS,
    output:
        variable_windows = "internal_output/DMN/homer/region_matched_background/variable/{libname}.{sample_label}.sampled_variable_windows.bed.gz",
        fixed_windows = "internal_output/DMN/homer/region_matched_background/fixed/{libname}.{sample_label}.sampled_fixed_windows.bed.gz"
    params:
        error_out_file = "stderr/{libname}.{sample_label}.sample_background_windows_by_region.err",
        out_file = "stdout/{libname}.{sample_label}.sample_background_windows_by_region.out",
        run_time = "10:00",
        memory = "3000",
        cores = 1,
        job_name = "sample_background_windows"
    benchmark: "benchmarks/sample_background_windows_by_region/{libname}.{sample_label}.sample_background_windows_by_region.txt"
    shell:
        """
        {R_EXE} --vanilla {SCRIPT_PATH}/sample_matched_background_by_region.R \
            {input.enriched_windows} \
            {input.all_windows} \
            75 \
            internal_output/DMN/homer/region_matched_background \
            {wildcards.libname}.{wildcards.sample_label};
        """

rule run_homer_dmn:
    input:
        finemapped_windows = "internal_output/DMN/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
        background = "internal_output/DMN/homer/region_matched_background/fixed/{libname}.{sample_label}.sampled_fixed_windows.bed.gz",
        genome = config['GENOMEFA']
    output:
        report = "internal_output/DMN/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
        motif = "internal_output/DMN/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerMotifs.all.motifs"
    params:
        error_out_file = "stderr/{libname}.{sample_label}.{signal_type}.run_homer.err",
        out_file = "stdout/{libname}.{sample_label}.{signal_type}.run_homer.out",
        run_time = "40:00",
        memory = "2000",
        cores = 1,
        job_name = "run_homer"
    conda:
         "envs/homer.yaml"
    benchmark: "benchmarks/run_homer/{libname}.{sample_label}.{signal_type}.txt"
    shell:
        """
        findMotifsGenome.pl <(less {input.finemapped_windows} | awk -v OFS=\"\t\" '{{print $4 \":\"$9,$1,$2+1,$3,$6}}') \
            {input.genome} internal_output/DMN/homer/finemapped_results/{wildcards.signal_type}/{wildcards.libname}.{wildcards.sample_label} \
            -preparsedDir internal_output/DMN/homer/preparsed -size given -rna -nofacts -S 20 -len 5,6,7,8,9 -nlen 1 \
            -bg <(zcat {input.background} | awk -v OFS=\"\t\" '{{print $4,$1,$2+1,$3,$6}}') 
        """

rule sample_background_windows_by_region_dmn_most_enriched:
    input:
        enriched_windows = "internal_output/DMN/most_enriched_selected/{libname}.{sample_label}.enriched_window.tsv",
        all_windows = FEATURE_ANNOTATIONS,
    output:
        variable_windows = "internal_output/DMN/homer_most_enriched/region_matched_background/variable/{libname}.{sample_label}.sampled_variable_windows.bed.gz",
        fixed_windows = "internal_output/DMN/homer_most_enriched/region_matched_background/fixed/{libname}.{sample_label}.sampled_fixed_windows.bed.gz"
    params:
        error_out_file = "stderr/{libname}.{sample_label}.sample_background_windows_by_region.err",
        out_file = "stdout/{libname}.{sample_label}.sample_background_windows_by_region.out",
        run_time = "10:00",
        memory = "3000",
        cores = 1,
        job_name = "sample_background_windows",
        outdir = "internal_output/DMN/homer_most_enriched"
    benchmark: "benchmarks/sample_background_windows_by_region/{libname}.{sample_label}.sample_background_windows_by_region.txt"
    shell:
        """
        {R_EXE} --vanilla {SCRIPT_PATH}/sample_matched_background_by_region.R \
            {input.enriched_windows} \
            {input.all_windows} \
            75 \
            {params.outdir}/region_matched_background \
            {wildcards.libname}.{wildcards.sample_label};
        """

rule run_homer_dmn_most_enriched:
    input:
        finemapped_windows = "internal_output/DMN/finemapping/mapped_sites_most_enriched/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
        background = "internal_output/DMN/homer_most_enriched/region_matched_background/fixed/{libname}.{sample_label}.sampled_fixed_windows.bed.gz",
        genome = config['GENOMEFA']
    output:
        report = "internal_output/DMN/homer_most_enriched/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
        motif = "internal_output/DMN/homer_most_enriched/finemapped_results/{signal_type}/{libname}.{sample_label}/homerMotifs.all.motifs"
    params:
        error_out_file = "stderr/{libname}.{sample_label}.{signal_type}.run_homer.err",
        out_file = "stdout/{libname}.{sample_label}.{signal_type}.run_homer.out",
        run_time = "40:00",
        memory = "2000",
        cores = 1,
        job_name = "run_homer",
        outdir = "internal_output/DMN/homer_most_enriched"
    conda:
         "envs/homer.yaml"
    benchmark: "benchmarks/run_homer/{libname}.{sample_label}.{signal_type}.txt"
    shell:
        """
        findMotifsGenome.pl <(less {input.finemapped_windows} | awk -v OFS=\"\t\" '{{print $4 \":\"$9,$1,$2+1,$3,$6}}') \
            {input.genome} {params.outdir}/finemapped_results/{wildcards.signal_type}/{wildcards.libname}.{wildcards.sample_label} \
            -preparsedDir {params.outdir}/preparsed -size given -rna -nofacts -S 20 -len 5,6,7,8,9 -nlen 1 \
            -bg <(zcat {input.background} | awk -v OFS=\"\t\" '{{print $4,$1,$2+1,$3,$6}}') 
        """

# rule run_homer:
#     input:
#         enriched_windows = 'homer/inputs/{libname}.{clip_sample_label}.{region}.enrich.homer',
#         background = 'homer/inputs/{libname}.{clip_sample_label}.{region}.tested.homer',
#         genome = config['GENOME_FA']
#     output:
#         report = "output/homer/results/{libname}.{clip_sample_label}.{region}/homerResults.html",
#         motif = "output/homer/results/{libname}.{clip_sample_label}.{region}/homerMotifs.all.motifs",
#     params:
#         error_out_file = "stderr/{libname}.{clip_sample_label}.{region}.run_homer.err",
#         out_file = "stdout/{libname}.{clip_sample_label}.{region}.run_homer.out",
#         run_time = "40:00",
#         memory = "2000",
#         job_name = "run_homer",
#         cores = 1
#     benchmark: "benchmarks/run_homer/{libname}.{clip_sample_label}.{region}.all_replicates.reproducible.txt"
#     conda:
#         "envs/homer.yaml"
#     shell:
#         """
#         n_enrich=$(cat {input.enriched_windows} | wc -l)
#         if [[ "$n_enrich" -le 50 ]]
#         then
#             echo "less than 50 region" > {output.motif}
#             touch {output.report}
#         else
#             findMotifsGenome.pl {input.enriched_windows} \
#             {input.genome} \
#             output/homer/results/{wildcards.libname}.{wildcards.clip_sample_label}.{wildcards.region} \
#             -bg {input.background} \
#             -preparsedDir output/homer/preparsed \
#             -size given \
#             -rna \
#             -nofacts \
#             -len 5,6,7,8,9 \
#             -nlen 1 
#         fi
#         """