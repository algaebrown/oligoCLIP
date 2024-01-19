rule annotate:
    # https://github.com/byee4/annotator
    input:
        peak="{something}.normed.compressed.bed"
    output:
        "{something}.normed.compressed.annotate.bed"
    params:
        run_time="3:00:00",
        error_out_file = "error_files/annotate.{something}",
        out_file =  "stdout/annotate.{something}",
        cores = "1",
        gtf_db = config['GTF'],
        sps = config['ANNOTATOR_SPECIES'],
        memory = 40000,
    shell:
        """
        module load annotator/0.0.15
        annotator --input {input.peak} \
        --output {output} \
        --gtfdb {params.gtf_db} \
        --species {params.sps}
        """
        
rule motif_analysis:
    # https://github.com/YeoLab/clip_analysis_legacy
    input:
        peak="{something}.normed.compressed.bed"
    output:
        filtered_peak = "{something}.normed.compressed.filtered.bed",
        pickle="{something}.normed.compressed.motif.pickle",
        svg="{something}.normed.compressed.motif.svg",
    params:
        run_time="3:00:00",
        error_out_file = "error_files/motif.{something}",
        out_file =  "stdout/motif.{something}",
        cores = "1",
        fa = config['GENOMEFA'],
        sps = config['ANNOTATOR_SPECIES'],
        homer="{something}.normed.compressed.homer",
        log10p_thres=3,
        l2fc_thres = 3,
        memory = 80000,
    shell:
        """
        module load eclipanalysis

        awk '{{ if ($4 > {params.log10p_thres})&($5 > {params.l2fc_thres}) print }}' {input.peak}  > {output.filtered_peak}

        analyze_motifs --peaks {input.peak} --k 6 \
            --out_pickle_file {output.pickle} \
            --out_homer_dir {params.homer} \
            --out_file {output.svg} \
            --species {params.sps} \
            --genome_fasta {params.fa}
        """