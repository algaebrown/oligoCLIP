rule extract_read_two:
    input:
        #bam="input/{sample_label}.bam"
        bam=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["bam"].values[0]),
    output:
        read2="processed_bam/{sample_label}.r2.bam",
        read1="processed_bam/{sample_label}.r1.bam"
    params:
        run_time="2:00:00",
        cores = 1,
        error_out_file = "error_files/extract_read2",
        out_file = "stdout/extract_read2.{sample_label}"
    shell:
        """
        module load samtools
        paired=$(samtools view -c -f 1 {input.bam})
        if [ "$paired" -lt "1" ]; 
            then cp {input.bam} {output.read2};
            cp {input.bam} {output.read1}
        else


            # # get read2 only
            samtools view -h -f 0x0080 {input.bam} | samtools view -Sb - > {output.read2}

            # get read1 only
            samtools view -h -f 0x0040 {input.bam} | samtools view -Sb - > {output.read1}
        fi
        """
rule CITS_bam_to_bedgragh:
    input:
        bam="processed_bam/{sample_label}.r2.bam"
    output:
        pos="CITS/{sample_label}.pos.bedgraph",
        neg="CITS/{sample_label}.neg.bedgraph"
    params:
        run_time="2:00:00",
        error_out_file = "error_files/CITS_bedgraph",
        cores = 1,
        out_file = "stdout/CITS_bedgraph.{sample_label}"
    shell:
        """
        module load bedtools;
        set +o pipefail;
        
        bedtools genomecov -ibam {input.bam} -bg -scale 1 -strand + -5 > {output.pos}
        bedSort {output.pos} {output.pos}

        bedtools genomecov -ibam {input.bam} -bg -scale -1 -strand - -5 > {output.neg}
        bedSort {output.neg} {output.neg}
        """
        
rule COV_bam_to_bedgraph:
    input:
        bam="processed_bam/{sample_label}.r2.bam"
    output:
        pos="coverage/{sample_label}.pos.bedgraph",
        neg="coverage/{sample_label}.neg.bedgraph"
    params:
        run_time="1:00:00",
        error_out_file = "error_files/coverage_bedgraph",
        out_file = "stdout/COV_bedgraph.{sample_label}",
        cores = 1,
    shell:
        """
        module load bedtools;
        set +o pipefail;

        # coverage
        bedtools genomecov -ibam {input.bam} -bg -scale 1 -strand + -split > {output.pos}
        bedSort {output.pos} {output.pos}
        bedtools genomecov -ibam {input.bam} -bg -scale -1 -strand - -split > {output.neg}
        bedSort {output.neg} {output.neg}
        """

rule bedgraph_to_bw:
    input:
        bedgraph="{something}.bedgraph",
    output:
        bw="{something}.bw"
    params:
        run_time="6:00:00",
        chr_size=config['CHROM_SIZES'],
        error_out_file = "error_files/{something}.bedgraph_to_bw",
        out_file = "stdout/{something}.bedgraph_to_bw",
        cores = 1,
    shell:
        """
        module load ucsctools
        bedGraphToBigWig {input.bedgraph} {params.chr_size} {output.bw}
        """