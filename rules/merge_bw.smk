rule merge_other_bw_as_bg:
    input:
        bws=lambda wildcards: expand("{libname}/bw/{signal_type}/{sample_label}.{strand}.bw",
            libname = [wildcards.libname],
            sample_label = list(set(config['rbps'])-set([wildcards.sample_label])-set(config['AS_INPUT'])),
            signal_type = [wildcards.signal_type],
            strand = [wildcards.strand]
        )
    output:
        tem=temp("{libname}/bw_bg/{signal_type}/{sample_label}.{strand}.temp.bedgraph"),
        merged_bedgraph=temp("{libname}/bw_bg/{signal_type}/{sample_label}.{strand}.bedgraph"),
    params:
        run_time="2:00:00",
        error_out_file = "error_files/merge_bw.{libname}.{signal_type}.{sample_label}.{strand}.err",
        out_file = "stdout/merge_bw.{libname}.{signal_type}.{sample_label}.{strand}.out",
        cores = 1,
        memory = 80000,
    conda:
        "envs/bwmerge.yaml"
    shell:
        """
        bigWigMerge {input.bws} {output.tem} -threshold=-99999999
        bedSort {output.tem} {output.merged_bedgraph}
        """