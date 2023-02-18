SCRIPT_PATH = config['SCRIPT_PATH']
R_EXE = config['R_EXE']
rule get_nt_coverage:
    input:
        windows = "internal_output/enriched_windows/{libname}.{sample_label}.enriched_windows.tsv.gz",
        clip_bw_pos = "{libname}/bw/{signal_type}/{sample_label}.pos.bw",
        clip_bw_neg = "{libname}/bw/{signal_type}/{sample_label}.neg.bw",
        input_bw_pos = "{libname}/bw_bg/{signal_type}/{sample_label}.pos.bw",
        input_bw_neg = "{libname}/bw_bg/{signal_type}/{sample_label}.neg.bw",
    output:
        nt_coverage = "internal_output/finemapping/nt_coverage/{signal_type}/{libname}.{sample_label}.nt_coverage.bed"
    params:
        error_out_file = "stderr/{libname}.{sample_label}.{signal_type}.get_nt_coverage.err",
        out_file = "stdout/{libname}.{sample_label}.{signal_type}.get_nt_coverage.out",
        run_time = "1:00:00",
        memory = "15000",
        job_name = "get_nt_coverage",
        cores = 1
    conda:
        "envs/metadensity.yaml"
    benchmark: "benchmarks/get_nt_coverage/{libname}.{sample_label}.{signal_type}.all_replicates.reproducible.txt"
    shell:
        """
        python {SCRIPT_PATH}/prepare_finemap.py \
        --ipminus {input.clip_bw_neg} \
        --ipplus {input.clip_bw_pos} \
        --inminus {input.input_bw_neg} \
        --inplus {input.input_bw_pos} \
        --region {input.windows} \
        --bed {output.nt_coverage}
        """
rule finemap_windows:
    input:
        nt_coverage = "internal_output/finemapping/nt_coverage/{signal_type}/{libname}.{sample_label}.nt_coverage.bed",
        # col_names = c("chr","start","end","name","score","strand","window_n","input","clip")        
    output:
        finemapped_windows = "internal_output/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz"
    params:
        error_out_file = "error_files/{libname}.{sample_label}.{signal_type}.finemap_windows.err",
        out_file = "stdout/{libname}.{sample_label}.{signal_type}.finemap_windows.out",
        run_time = "1:00:00",
        memory = "10000",
        job_name = "finemap_windows",
        cores = 1,
        outdir = "internal_output/finemapping/mapped_sites/{signal_type}"
    benchmark: "benchmarks/finemap_windows/{libname}.{sample_label}.{signal_type}.all_replicates.reproducible.txt"
    shell:
        """
        {R_EXE} --vanilla {SCRIPT_PATH}/finemap_enriched_windows.R \
            {input.nt_coverage} \
            {params.outdir} \
            {wildcards.libname}.{wildcards.sample_label}
        """

####### DMN windows ########
rule get_nt_coverage_dmn:
    input:
        windows = "internal_output/DMN/{libname}.{sample_label}.enriched_window.tsv",
        clip_bw_pos = "{libname}/bw/{signal_type}/{sample_label}.pos.bw",
        clip_bw_neg = "{libname}/bw/{signal_type}/{sample_label}.neg.bw",
        input_bw_pos = "{libname}/bw_bg/{signal_type}/{sample_label}.pos.bw",
        input_bw_neg = "{libname}/bw_bg/{signal_type}/{sample_label}.neg.bw",
    output:
        nt_coverage = "internal_output/DMN/finemapping/nt_coverage/{signal_type}/{libname}.{sample_label}.nt_coverage.bed"
    params:
        error_out_file = "stderr/{libname}.{sample_label}.{signal_type}.get_nt_coverage.err",
        out_file = "stdout/{libname}.{sample_label}.{signal_type}.get_nt_coverage.out",
        run_time = "1:00:00",
        memory = "15000",
        job_name = "get_nt_coverage",
        cores = 1
    conda:
        "envs/metadensity.yaml"
    benchmark: "benchmarks/DMN/get_nt_coverage.{libname}.{sample_label}.{signal_type}"
    shell:
        """
        python {SCRIPT_PATH}/prepare_finemap.py \
        --ipminus {input.clip_bw_neg} \
        --ipplus {input.clip_bw_pos} \
        --inminus {input.input_bw_neg} \
        --inplus {input.input_bw_pos} \
        --region {input.windows} \
        --bed {output.nt_coverage}
        """

rule finemap_windows_dmn:
    input:
        nt_coverage = "internal_output/DMN/finemapping/nt_coverage/{signal_type}/{libname}.{sample_label}.nt_coverage.bed",
    output:
        finemapped_windows = "internal_output/DMN/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz"
    params:
        error_out_file = "error_files/{libname}.{sample_label}.{signal_type}.finemap_windows.err",
        out_file = "stdout/{libname}.{sample_label}.{signal_type}.finemap_windows.out",
        run_time = "1:00:00",
        memory = "10000",
        job_name = "finemap_windows",
        cores = 1,
        outdir = "internal_output/DMN/finemapping/mapped_sites/{signal_type}"
    benchmark: "benchmarks/DMN/finemap.{libname}.{sample_label}.{signal_type}"
    shell:
        """
        {R_EXE} --vanilla {SCRIPT_PATH}/finemap_enriched_windows.R \
            {input.nt_coverage} \
            {params.outdir} \
            {wildcards.libname}.{wildcards.sample_label}
        """

rule get_nt_coverage_dmn_most_enriched:
    input:
        windows = "internal_output/DMN/most_enriched_selected/{libname}.{sample_label}.enriched_window.tsv",
        clip_bw_pos = "{libname}/bw/{signal_type}/{sample_label}.pos.bw",
        clip_bw_neg = "{libname}/bw/{signal_type}/{sample_label}.neg.bw",
        input_bw_pos = "{libname}/bw_bg/{signal_type}/{sample_label}.pos.bw",
        input_bw_neg = "{libname}/bw_bg/{signal_type}/{sample_label}.neg.bw",
    output:
        nt_coverage = "internal_output/DMN/finemapping/nt_coverage_most_enriched/{signal_type}/{libname}.{sample_label}.nt_coverage.bed"
    params:
        error_out_file = "stderr/{libname}.{sample_label}.{signal_type}.get_nt_coverage.err",
        out_file = "stdout/{libname}.{sample_label}.{signal_type}.get_nt_coverage.out",
        run_time = "1:00:00",
        memory = "15000",
        job_name = "get_nt_coverage",
        cores = 1
    conda:
        "envs/metadensity.yaml"
    benchmark: "benchmarks/DMN/get_nt_coverage.{libname}.{sample_label}.{signal_type}"
    shell:
        """
        python {SCRIPT_PATH}/prepare_finemap.py \
        --ipminus {input.clip_bw_neg} \
        --ipplus {input.clip_bw_pos} \
        --inminus {input.input_bw_neg} \
        --inplus {input.input_bw_pos} \
        --region {input.windows} \
        --bed {output.nt_coverage}
        """

rule finemap_windows_dmn_most_enriched:
    input:
        nt_coverage = "internal_output/DMN/finemapping/nt_coverage_most_enriched/{signal_type}/{libname}.{sample_label}.nt_coverage.bed"
    output:
        finemapped_windows = "internal_output/DMN/finemapping/mapped_sites_most_enriched/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz"
    params:
        error_out_file = "error_files/{libname}.{sample_label}.{signal_type}.finemap_windows.err",
        out_file = "stdout/{libname}.{sample_label}.{signal_type}.finemap_windows.out",
        run_time = "1:00:00",
        memory = "10000",
        job_name = "finemap_windows",
        cores = 1,
        outdir = "internal_output/DMN/finemapping/mapped_sites_most_enriched/{signal_type}"
    benchmark: "benchmarks/DMN/finemap.{libname}.{sample_label}.{signal_type}"
    shell:
        """
        {R_EXE} --vanilla {SCRIPT_PATH}/finemap_enriched_windows.R \
            {input.nt_coverage} \
            {params.outdir} \
            {wildcards.libname}.{wildcards.sample_label}
        """

rule finemap_to_bedgraph:
    input:
        "{something}.bed.gz"
    output:
        pos = "{something}.pos.bedgraph",
        neg = "{something}.neg.bedgraph"
    params:
        error_out_file = "error_files/bed2bedgraph.{something}.err",
        out_file = "stdout/bed2bedgraph.{something}.out",
        run_time = "00:10:00",
        memory = "10000",
        job_name = "peak_to_bedgraph",
        cores = 1,
    shell:
        """
        zcat {input} | grep "-" | awk '{{ print $1"\t"$2"\t"$3"\t"1 }}' > {output.neg}
        zcat {input} | grep "+" | awk '{{ print $1"\t"$2"\t"$3"\t"1 }}' > {output.pos}
        """