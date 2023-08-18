from pathlib import Path
SCRIPT_PATH = config['SCRIPT_PATH']
R_EXE = config['R_EXE']


rule get_nt_coverage:
    input:
        windows = lambda wildcards: f"{wildcards.root_dir}/enriched_windows/{wildcards.libname}.{wildcards.sample_label}.enriched_windows.tsv.gz" if 'skipper' 
            in wildcards.root_dir else f"{wildcards.root_dir}/{wildcards.libname}.{wildcards.sample_label}.enriched_windows.tsv",
        clip_bw_pos = "{libname}/bw/{signal_type}/{sample_label}.pos.bw",
        clip_bw_neg = "{libname}/bw/{signal_type}/{sample_label}.neg.bw",
        input_bw_pos = lambda wildcards: f"{wildcards.libname}/bw_bg/{wildcards.signal_type}/{wildcards.sample_label}.pos.bw" if 'external' not in wildcards.root_dir 
            else f"external_bw/{wildcards.signal_type}/"+wildcards.root_dir.split('/')[-1]+".pos.bw",
        input_bw_neg = lambda wildcards: f"{wildcards.libname}/bw_bg/{wildcards.signal_type}/{wildcards.sample_label}.neg.bw" if 'external' not in wildcards.root_dir 
            else f"external_bw/{wildcards.signal_type}/"+wildcards.root_dir.split('/')[-1]+".neg.bw"
    output:
        nt_coverage = "{root_dir}/finemapping/nt_coverage/{signal_type}/{libname}.{sample_label}.nt_coverage.bed"
    params:
        error_out_file = "error_files/get_nt_coverage.{root_dir}.{libname}.{sample_label}.{signal_type}.err",
        out_file = "stdout/get_nt_coverage.{root_dir}.{libname}.{sample_label}.{signal_type}.out",
        run_time = "1:00:00",
        memory = "15000",
        job_name = "get_nt_coverage",
        cores = 1
    conda:
        "envs/metadensity.yaml"
    benchmark: "benchmarks/get_nt_coverage/{root_dir}.{libname}.{sample_label}.{signal_type}.txt"
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
        nt_coverage = "{root_dir}/finemapping/nt_coverage/{signal_type}/{libname}.{sample_label}.nt_coverage.bed",
        # col_names = c("chr","start","end","name","score","strand","window_n","input","clip")        
    output:
        finemapped_windows = "{root_dir}/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz"
    params:
        error_out_file = "error_files/finemap_windows.{root_dir}.{libname}.{sample_label}.{signal_type}.err",
        out_file = "stdout/finemap_windows.{root_dir}.{libname}.{sample_label}.{signal_type}.out",
        run_time = "1:00:00",
        memory = "10000",
        job_name = "finemap_windows",
        cores = 1,
        outdir = lambda wildcards, output: Path(output.finemapped_windows).parent
    benchmark: "benchmarks/finemap_windows/{root_dir}.{libname}.{sample_label}.{signal_type}.all_replicates.reproducible.txt"
    shell:
        """
        {R_EXE} --vanilla {SCRIPT_PATH}/finemap_enriched_windows.R \
            {input.nt_coverage} \
            {params.outdir} \
            {wildcards.libname}.{wildcards.sample_label}
        """

# ####### beta-mixture/DMM windows ########
# use rule get_nt_coverage as get_nt_coverage_beta_mixture with:
#     input:
#         windows = "{root_dir}/{libname}.{sample_label}.enriched_windows.tsv",
#         clip_bw_pos = "{libname}/bw/{signal_type}/{sample_label}.pos.bw",
#         clip_bw_neg = "{libname}/bw/{signal_type}/{sample_label}.neg.bw",
#         input_bw_pos = "{libname}/bw_bg/{signal_type}/{sample_label}.pos.bw",
#         input_bw_neg = "{libname}/bw_bg/{signal_type}/{sample_label}.neg.bw",
#     output:
#         nt_coverage = "{root_dir}/finemapping/nt_coverage/{signal_type}/{libname}.{sample_label}.nt_coverage.bed"

# use rule finemap_windows as finemap_windows_beta with:
#     input:
#         nt_coverage = "{mixture_root_dir}/finemapping/nt_coverage/{signal_type}/{libname}.{sample_label}.nt_coverage.bed",
#     output:
#         finemapped_windows = "{mixture_root_dir}/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz"
    
# ########### DMM #############
# use rule get_nt_coverage as get_nt_coverage_dmm with:
#     input:
#         windows = "DMM/{libname}.{sample_label}.enriched_windows.tsv",
#         clip_bw_pos = "{libname}/bw/{signal_type}/{sample_label}.pos.bw",
#         clip_bw_neg = "{libname}/bw/{signal_type}/{sample_label}.neg.bw",
#         input_bw_pos = "{libname}/bw_bg/{signal_type}/{sample_label}.pos.bw",
#         input_bw_neg = "{libname}/bw_bg/{signal_type}/{sample_label}.neg.bw",
#     output:
#         nt_coverage = "DMM/finemapping/nt_coverage/{signal_type}/{libname}.{sample_label}.nt_coverage.bed"

# use rule finemap_windows as finemap_windows_dmm with:
#     input:
#         nt_coverage = "DMM/finemapping/nt_coverage/{signal_type}/{libname}.{sample_label}.nt_coverage.bed",
#     output:
#         finemapped_windows = "DMM/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz"

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