import pandas as pd
locals().update(config)
rbps = config['rbps']
experiments = config['experiments']
libnames = config['libnames']

manifest = pd.read_table(config['MANIFEST'], index_col = False, sep = ',')

def experiment_to_libname(experiment):
    libnames = manifest.loc[manifest['experiment']==experiment, 'libname'].tolist()
    assert len(libnames)>0
    return libnames
def libname_to_experiment(libname):
    return manifest.loc[manifest['libname']==libname, 'experiment'].iloc[0]


rule uniq_repeats:
    input:
        repeatmasker = config['REPEAT_TABLE'],
        genome = config['GENOMEFA']
    output:
        sorted_bed = temp("repeats.sort.temp.bed.gz"),
        unique_repeats = config['REPEAT_TABLE'].replace(".tsv", ".sort.unique.bed")
    params:
        error_out_file = "error_files/calc_partition_nuc.err",
        out_file = "stdout/calc_partition_nuc.out",
        run_time = "40:00",
        memory = 1000,
        cores = 1,
    benchmark: "benchmarks/uniq_repeats.txt"
    conda:
        "envs/bedtools.yaml"
    shell:
        "zcat {REPEAT_TABLE} | awk -v OFS=\"\\t\" '{{print $6,$7,$8,$11 \":\" name_count[$11]++, $2, $10,$11,$12,$13}} "
            "$13 == \"L1\" || $13 == \"Alu\" {{$11 = $11 \"_AS\"; $12 = $12 \"_AS\"; $13 = $13 \"_AS\"; "
            "if($10 == \"+\") {{$10 = \"-\"}} else {{$10 = \"+\"}}; print $6,$7,$8,$11 \":\" name_count[$11]++, $2, $10,$11,$12,$13}}' | "
            "tail -n +2 | bedtools sort -i - | gzip > {output.sorted_bed}; "
        "bedtools coverage -s -d -a {output.sorted_bed} -b {output.sorted_bed}  | awk -v OFS=\"\\t\" "
            "'$NF >1 {{print $1,$2+$(NF-1)-1,$2+$(NF-1),$4,$5,$6}}' | "
            "bedtools sort -i - | "
            "bedtools merge -c 4,5,6 -o distinct -s -i - | "
            "bedtools subtract -s -a {output.sorted_bed} -b - | "
            "bedtools nuc -s -fi {input.genome} -bed -  | awk -v OFS=\"\\t\" 'NR > 1 {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11}}' | "
            "gzip -c > {output.unique_repeats}"

rule quantify_repeats:
    input:
        config['CHROM_SIZES'],
        bam = "{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",
        repeats = rules.uniq_repeats.output.unique_repeats,
        #repeats = config['REPEAT_TABLE'].replace(".tsv", ".sort.unique.bed")
    output:
        counts = "counts/repeats/vectors/{libname}.{sample_label}.counts"
    params:
        error_out_file = "error_files/quantify_repeats.{libname}.{sample_label}.err",
        out_file = "stdout/quantify_repeats.{libname}.{sample_label}.out",
        run_time = "15:00",
        memory = 20000,
        replicate_label = "{libname}.{sample_label}",
        cores = 1,
        uninformative_read = config['UNINFORMATIVE_READ']
    benchmark: "benchmarks/repeats/unassigned_experiment.{libname}.{sample_label}.quantify_repeats.txt"
    # container:
    #     "docker://howardxu520/skipper:bigwig_1.0" # keep having problems with TSCC 2.0 with docker
    conda:
        "envs/bedtools.yaml"
    shell:
        "bedtools bamtobed -i {input.bam} | awk '($1 != \"chrEBV\") && ($4 !~ \"/{params.uninformative_read}$\")' | "
            "bedtools flank -s -l 1 -r 0 -g {CHROM_SIZES} -i - | "
            "bedtools shift -p 1 -m -1 -g {CHROM_SIZES} -i - | "
            "bedtools sort -i - | "
            "bedtools coverage -s -counts -a {input.repeats} -b - | "
            "awk 'BEGIN {{print \"{params.replicate_label}\"}} {{print $NF}}' > {output.counts}"


rule count_repeat_tables:
    input:
        unique_repeats = config['REPEAT_TABLE'].replace(".tsv", ".sort.unique.bed"),
        replicate_counts = lambda wildcards: expand(
            "counts/repeats/vectors/{libname}.{sample_label}.counts", 
            libname = experiment_to_libname(wildcards.experiment), # TODO: make dictionary
            sample_label = [wildcards.sample_label]),
    output:
        name_table = "counts/repeats/tables/name/{experiment}.{sample_label}.tsv.gz",
        class_table = "counts/repeats/tables/class/{experiment}.{sample_label}.tsv.gz",
        family_table = "counts/repeats/tables/family/{experiment}.{sample_label}.tsv.gz",
    params:
        error_out_file = "error_files/make_repeat_count_tables.{experiment}.{sample_label}.err",
        out_file = "stdout/make_repeat_count_tables.{experiment}.{sample_label}..out",
        run_time = "00:15:00",
        cores = "1",
        memory = 200,
    benchmark: "benchmarks/counts/{experiment}.{sample_label}.all_replicates.make_repeat_count_table.txt"
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

rule fit_clip_betabinomial_re_model:
    input:
        table = lambda wildcards: "counts/repeats/tables/name/"+libname_to_experiment(wildcards.libname)+".{sample_label}.tsv.gz",
    output:
        coef = "skipper_CC/clip_model_coef_re/{libname}.{sample_label}.tsv",
        # plot = lambda wildcards: expand("output/figures/clip_distributions/{{experiment_label}}.{{clip_replicate_label}}.{other_label}.clip_distribution.pdf", other_label = experiment_to_input_replicate_labels[wildcards.experiment_label][wildcards.Input_replicate_label])
    params:
        error_out_file = "error_files/fit_clip_betabinomial_re_model.{libname}.{sample_label}.err",
        out_file = "stdout/fit_clip_betabinomial_re_model.{libname}.{sample_label}.out",
        run_time = "10:00",
        memory = 40000,
        cores = 1,
    benchmark: "benchmarks/fit_clip_betabinomial_re_model/{libname}.{sample_label}.fit_clip.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/fit_clip_betabinom_re.R \
            {input.table} \
            {wildcards.libname} \
            {wildcards.libname}.{wildcards.sample_label} \
            {output.coef}
        """

rule sum_all_other_background_re:
    input:
        lambda wildcards: expand("counts/repeats/vectors/{libname}.{sample_label}.counts",
        libname = ["{libname}"],
        sample_label = list(set(rbps)-set([wildcards.sample_label])-set([config['AS_INPUT']]))
        )
    output:
        counts= "counts_CC/repeats/vectors/{libname}.{sample_label}.counts",
    params:
        error_out_file = "error_files/sum_re_reads.{libname}.{sample_label}.err",
        out_file = "stdout/sum_re_reads.{libname}.{sample_label}.out",
        run_time = "20:00",
        cores = "1",
        memory = 10000,
        replicate_label = "{libname}.internal",
    benchmark: "benchmarks/counts/unassigned_experiment.{libname}.{sample_label}.sum_read.txt"
    shell:
        """
        awk '{{arr[FNR]+=$1}}END{{for(i=2;i<=FNR;i+=1){{print arr[i]}} }}' {input} | awk 'BEGIN {{print \"{params.replicate_label}\"}} {{print}}' > {output}
        """

rule combine_ip_to_internal_background: # bg_sample_label = 'internal'
    input:
        count_table = "counts/repeats/tables/name/{experiment}.{sample_label}.tsv.gz",
        bg_counts = lambda wildcards: expand("counts_CC/repeats/vectors/{libname}.{sample_label}.counts", 
            libname = experiment_to_libname(wildcards.experiment), # TODO: make dictionary
            sample_label = [wildcards.sample_label])
    output:
        combined_count_table = "counts_CC/repeats/bgtables/internal/{experiment}.{sample_label}.tsv.gz"
    params:
        error_out_file = "error_files/combine_ip_to_background.{experiment}.{sample_label}.combine.err",
        out_file = "stdout/combine_ip_to_background.{experiment}.{sample_label}.combine.out",
        run_time = "1:00:00",
        cores = "1",
        memory = 10000,
    benchmark: "benchmarks/combine_table/{experiment}.internal.{sample_label}.combine.txt"
    shell:
        """
        paste <(zcat {input.count_table}) {input.bg_counts}| gzip -c > {output.combined_count_table}
        """

rule call_enriched_re:
    input:
        table = lambda wildcards: f"counts_CC/repeats/bgtables/internal/"+libname_to_experiment(wildcards.libname)+f".{wildcards.sample_label}.tsv.gz",
        repeats = config['REPEAT_TABLE'].replace(".tsv", ".sort.unique.bed"),
        parameters = "skipper_CC/clip_model_coef_re/{libname}.{sample_label}.tsv",
    output:
        "skipper_CC/figures/clip_scatter_re/{libname}.{sample_label}.clip_test_distribution.pdf",
        "skipper_CC/enriched_re/{libname}.{sample_label}.enriched_re.tsv.gz"
    params:
        input_replicate_label = lambda wildcards: wildcards.libname+'.internal',
        error_out_file = "error_files/{libname}.{sample_label}.call_enriched_re.err",
        out_file = "stdout/{libname}.{sample_label}.call_enriched_re.out",
        run_time = "00:25:00",
        memory = 3000,
        cores = 1,
    benchmark: "benchmarks/call_enriched_re/{libname}.{sample_label}.call_enriched_re.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/call_enriched_re.R \
            {input.table} \
            {input.repeats} \
            {input.parameters} \
            {params.input_replicate_label} \
            {wildcards.libname}.{wildcards.sample_label} \
            {wildcards.libname}.{wildcards.sample_label} \
            skipper_CC
        """

# rule find_reproducible_enriched_re:
#     input:
#         windows = lambda wildcards: expand("output/enriched_re/{{experiment_label}}.{clip_replicate_label}.enriched_re.tsv.gz", clip_replicate_label = experiment_to_clip_replicate_labels[wildcards.experiment_label])
#     output:
#         reproducible_windows = "output/reproducible_enriched_re/{experiment_label}.reproducible_enriched_re.tsv.gz",
#     params:
#         error_file = "stderr/{experiment_label}.find_reproducible_enriched_re.err",
#         out_file = "stdout/{experiment_label}.find_reproducible_enriched_re.out",
#         run_time = "5:00",
#         memory = "1000",
#     benchmark: "benchmarks/find_reproducible_enriched_re/{experiment_label}.all_replicates.reproducible.txt"
#     shell:
#         "Rscript --vanilla {TOOL_DIR}/identify_reproducible_re.R output/enriched_re/ {wildcards.experiment_label}"