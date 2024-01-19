locals().update(config)
rule make_repeat_megatable:
    input:
        unique_repeats = config['REPEAT_TABLE'].replace(".tsv", ".sort.unique.bed"),
        replicate_counts = lambda wildcards: expand(
            "counts/repeats/vectors/{libname}.{sample_label}.counts", 
            libname = ['{libname}'], # TODO: make dictionary
            sample_label = config['rbps']),
    output:
        name_table = "counts/repeats/megatables/name/{libname}.tsv.gz",
        class_table = "counts/repeats/megatables/class/{libname}.tsv.gz",
        family_table = "counts/repeats/megatables/family/{libname}.tsv.gz",
    params:
        error_out_file = "error_files/make_repeat_megatable.{libname}.err",
        out_file = "stdout/make_repeat_megatable.{libname}.out",
        run_time = "00:45:00",
        cores = "1",
        memory = 80000,
    benchmark: "benchmarks/counts/make_repeat_megatable.{libname}.txt"
    container: None
    shell:
        "echo \"repeat_name\" | paste - {input.replicate_counts} | sed -n '1p' | gzip > {output.name_table};"
        "echo \"repeat_class\" | paste - {input.replicate_counts} | sed -n '1p' | gzip > {output.class_table};"
        "echo \"repeat_family\" | paste - {input.replicate_counts} | sed -n '1p' | gzip > {output.family_table};"
        "paste <(zcat {input.unique_repeats} | awk -v OFS=\"\\t\" 'BEGIN {{print \"repeat_name\";}} {{print $7}}') {input.replicate_counts} | "
            "awk -v OFS=\"\\t\" 'NR > 1 {{for(i = 2; i <= NF; i++) {{tabulation[$1][i] += $i}} }} END {{for(name in tabulation) {{ printf name; for(i = 2; i <= NF; i++) {{printf \"\\t\" tabulation[name][i]}} print \"\";}} }}' | sort -k 1,1 | gzip >> {output.name_table};"
        "paste <(zcat {input.unique_repeats} | awk -v OFS=\"\\t\" 'BEGIN {{print \"repeat_class\";}} {{print $8}}') {input.replicate_counts} | "
            "awk -v OFS=\"\\t\" 'NR > 1 {{for(i = 2; i <= NF; i++) {{tabulation[$1][i] += $i}} }} END {{for(name in tabulation) {{ printf name; for(i = 2; i <= NF; i++) {{printf \"\\t\" tabulation[name][i]}} print \"\";}} }}' | sort -k 1,1 | gzip >> {output.class_table};"
        "paste <(zcat {input.unique_repeats} | awk -v OFS=\"\\t\" 'BEGIN {{print \"repeat_family\";}} {{print $9}}') {input.replicate_counts} | "
            "awk -v OFS=\"\\t\" 'NR > 1 {{for(i = 2; i <= NF; i++) {{tabulation[$1][i] += $i}} }} END {{for(name in tabulation) {{ printf name; for(i = 2; i <= NF; i++) {{printf \"\\t\" tabulation[name][i]}} print \"\";}} }}' | sort -k 1,1 | gzip >> {output.family_table};"

rule fit_DMN:
    input:
        table = "counts/repeats/megatables/{repeat_type}/{libname}.tsv.gz",
    output:
        "DMM_repeat/{repeat_type}/{libname}.goodness_of_fit.pdf",
        "DMM_repeat/{repeat_type}/{libname}.alpha.tsv",
        "DMM_repeat/{repeat_type}/{libname}.null.alpha.tsv",
        "DMM_repeat/{repeat_type}/{libname}.mixture_weight.tsv",
        "DMM_repeat/{repeat_type}/{libname}.weights.tsv"
    params:
        error_out_file = "error_files/fit_DMM_repeat.{libname}.{repeat_type}.err",
        out_file = "stdout/fit_DMM_repeat.{libname}.{repeat_type}.out",
        run_time = "12:00:00",
        memory = 10000,
        cores = "4",
        root_folder = lambda wildcards, output: Path(output[0]).parent
    benchmark: "benchmarks/DMM/fit.{libname}.{repeat_type}"
    conda:
        "envs/DMN.yaml"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/fit_DMN_multidimen_repeat.R \
            {input.table} \
            {params.root_folder} \
            {wildcards.libname} 
        """

rule analyze_DMN:
    input:
        "DMM_repeat/{repeat_type}/{libname}.goodness_of_fit.pdf",
        "DMM_repeat/{repeat_type}/{libname}.alpha.tsv",
        "DMM_repeat/{repeat_type}/{libname}.null.alpha.tsv",
        "DMM_repeat/{repeat_type}/{libname}.mixture_weight.tsv",
        "DMM_repeat/{repeat_type}/{libname}.weights.tsv",
        annotation = config['REPEAT_TABLE'],
        raw_counts = "counts/repeats/megatables/{repeat_type}/{libname}.tsv.gz",
        repeat_mask = 'mask/{libname}.repeat_mask.csv', #True means zscore > 1
    output:
        expand("DMM_repeat/{repeat_type}/{libname}.{sample_label}.enriched_windows.tsv", 
            sample_label = config['rbps'], 
            repeat_type = ['{repeat_type}'],
            libname = ['{libname}']),
        "DMM_repeat/{repeat_type}/{libname}.megaoutputs.tsv"
    params:
        error_out_file = "error_files/analyze_DMM_repeat.{libname}.{repeat_type}.err",
        out_file = "stdout/fit_DMM_repeat.{libname}.{repeat_type}.out",
        run_time = "1:00:00",
        memory = 10000,
        cores = "1",
    benchmark: "benchmarks/DMM/fit.{libname}.{repeat_type}"
    conda:
        "envs/tensorflow.yaml"
    shell:
        """
        python {SCRIPT_PATH}/analyze_DMM_repeat.py \
            {wildcards.libname} \
            . \
            {input.annotation}
        """
    
