rule merge_other_bw_as_bg:
    input:
        bws=lambda wildcards: expand("{libname}/bw/{signal_type}/{sample_label}.r1.{strand}.bw",
            libname = [wildcards.libname],
            sample_label = list(set(config['rbps'])-set([wildcards.sample_label])-set([config['AS_INPUT']])),
            signal_type = [wildcards.signal_type],
            strand = [wildcards.strand]
        )
    output:
        tem=temp("{libname}/bw_bg/{signal_type}/{sample_label}.r1.{strand}.temp.bedgraph"),
        merged_bedgraph="{libname}/bw_bg/{signal_type}/{sample_label}.r1.{strand}.bedgraph",
    params:
        run_time="2:00:00",
        error_out_file = "error_files/merge_bw.{libname}.{signal_type}.{sample_label}.{strand}",
        out_file = "stdout/{libname}.{signal_type}.{sample_label}.{strand}.bedgraph_to_bw",
        cores = 1,
    shell:
        """
        module load ucsctools
        bigWigMerge {input.bws} {output.tem} -threshold=-99999999
        bedSort {output.tem} {output.merged_bedgraph}
        """