from pathlib import Path
locals().update(config)

rule get_nt_coverage:
    input:
        windows = lambda wildcards: str(Path(wildcards.root_dir)/'enriched_windows'/+f"{wildcards.libname}.{wildcards.sample_label}.enriched_windows.tsv.gz") if 'skipper' 
            in wildcards.root_dir else Path(wildcards.root_dir)/(wildcards.libname+'.'+wildcards.sample_label+".enriched_windows.tsv"),
        clip_bw_pos = "{libname}/bw/{signal_type}/{sample_label}.pos.bw",
        clip_bw_neg = "{libname}/bw/{signal_type}/{sample_label}.neg.bw",
        input_bw_pos = lambda wildcards:Path(wildcards.libname)/'bw_bg'/wildcards.signal_type / (wildcards.sample_label+".pos.bw") if 'external' not in wildcards.root_dir 
            else f"external_bw/{wildcards.signal_type}/"+wildcards.root_dir.split('/')[-1]+".pos.bw",
        input_bw_neg = lambda wildcards:Path(wildcards.libname)/'bw_bg'/wildcards.signal_type / (wildcards.sample_label+".neg.bw") if 'external' not in wildcards.root_dir 
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
        outdir = lambda wildcards, output: str(Path(output.finemapped_windows).parent)
    benchmark: "benchmarks/finemap_windows/{root_dir}.{libname}.{sample_label}.{signal_type}.all_replicates.reproducible.txt"
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/finemap_enriched_windows.R \
            {input.nt_coverage} \
            {params.outdir} \
            {wildcards.libname}.{wildcards.sample_label}
        """

rule finemap_to_bedgraph:
    input:
        "{something}.bed.gz"
    output:
        pos = temp("{something}.pos.bedgraph"),
        neg = temp("{something}.neg.bedgraph")
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