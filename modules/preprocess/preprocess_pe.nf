#!/usr/bin/env nextflow

params.tile_length = 10
params.adapter_fwd = "AGATCGGAAGAGCACACGTC"
params.adapter_rev = "AGATCGGAAGAGCGTCGTGT"
params.main_path = "/projects/ps-yeolab3/bay001/codebase/oligoCLIP/"
params.quality_cutoff = 6
params.libName = 'testing'
params.umi_length = 10
params.outDir = '/projects/ps-yeolab3/bay001/codebase/oligoCLIP/nextflow_outputs'
params.greeting = 'Preprocessing eCLIP data!' 
params.barcode_csv = '/home/bay001/projects/codebase/oligoCLIP/test_files/iter7.csv'
params.reads = "/home/bay001/projects/codebase/oligoCLIP/test_files/GN_1020_SMALL.R{1,2}.fastq.gz"
params.star_sjdb = "/projects/ps-yeolab3/bay001/annotations/GRCh38/star_2_7_gencode40_sjdb/"
reads_ch = Channel.fromFilePairs(params.reads) 

greeting_ch = Channel.of(params.greeting) 

process tile_adapter {
    publishDir "${params.outDir}", mode: 'copy', overwrite: false
    input:
        val adapter_fwd
        val adapter_rev
    output:
        path "fwd_adapters.txt", emit: adapter_fwd_txt 
        path "rev_adapters.txt", emit: adapter_rev_txt
    script:
        """
        python ${params.main_path}/scripts/create_adaptor_tile.py ${params.adapter_fwd} fwd_adapters.txt ${params.tile_length}
        python ${params.main_path}/scripts/create_adaptor_tile.py ${params.adapter_rev} rev_adapters.txt ${params.tile_length}
        """
}
process trim_adapter {
    publishDir "${params.outDir}", mode: 'copy', overwrite: false
    input:
        val libName
        val adapter_fwd_txt
        val adapter_rev_txt
        tuple val(replicateId), path(reads)
    output:
        path "${libName}.Tr.R1.fastq.gz", emit: read1
        path "${libName}.Tr.R2.fastq.gz", emit: read2
        path "${libName}.Tr.metrics"
    script:
        """
        echo ${libName}
        echo ${adapter_fwd_txt}
        echo ${adapter_rev_txt}
        
        mkdir ${libName};
        cutadapt -a file:${adapter_fwd_txt} \
        -A file:${adapter_rev_txt} \
        --times 2 \
        -e 0.1 \
        --quality-cutoff ${params.quality_cutoff} \
        -m 23 \
        -o ${libName}.Tr.R1.fastq.gz \
        -p ${libName}.Tr.R2.fastq.gz \
        --cores=0 \
        ${reads[0]} ${reads[1]} > ${libName}.Tr.metrics
        """
}
process extract_umi_and_trim_polyG {
    publishDir "${params.outDir}", mode: 'copy', overwrite: false
    input:
        val libName
        path "${libName}.Tr.R1.fastq.gz"
        path "${libName}.Tr.R2.fastq.gz"
    output:
        path "${libName}.Tr.umi.R1.fastq.gz", emit: read1
        path "${libName}.Tr.umi.R2.fastq.gz", emit: read2
        path "${libName}.QC.umi.json"
        path "${libName}.QC.umi.html"

    script:
    """
    mkdir QC;
    fastp -i ${libName}.Tr.R1.fastq.gz \
        -I ${libName}.Tr.R2.fastq.gz \
        -o ${libName}.Tr.umi.R1.fastq.gz \
        -O ${libName}.Tr.umi.R2.fastq.gz \
        --disable_adapter_trimming \
        --umi \
        --umi_len=${params.umi_length} \
        --umi_loc=read1 \
        --trim_poly_g \
        -j ${libName}.QC.umi.json \
        -h ${libName}.QC.umi.html \
        -w ${task.cpus}
    """
}
process trim_umi_from_read2 {
    publishDir "${params.outDir}", mode: 'copy', overwrite: false
    input:
        val libName
        path "${libName}.Tr.umi.R1.fastq.gz"
        path "${libName}.Tr.umi.R2.fastq.gz"
    output:
        path("${libName}.Tr.umi.Tr.R1.fastq.gz"), emit: read1
        path("${libName}.Tr.umi.Tr.R2.fastq.gz"), emit: read2
    script:
    """
    mv ${libName}.Tr.umi.R1.fastq.gz ${libName}.Tr.umi.Tr.R1.fastq.gz
    zcat ${libName}.Tr.umi.R2.fastq.gz > ${libName}.Tr.umi.R2.fastq
    seqtk trimfq -e {params.umi_length} ${libName}.Tr.umi.R2.fastq | gzip > ${libName}.Tr.umi.Tr.R2.fastq.gz
    """
}
process demultiplex {
    publishDir "${params.outDir}", mode: 'copy', overwrite: false
    input:
        val libName
        path("${libName}.Tr.umi.R1.fastq.gz")
        path("${libName}.Tr.umi.R2.fastq.gz")
    output:
        path("ultraplex_demux_*.fastq.gz")
    script:
    """
    ultraplex -i ${libName}.Tr.umi.R2.fastq.gz -i2 ${libName}.Tr.umi.R1.fastq.gz -b ${params.barcode_csv}  \
        -m5 1 -m3 0 -t ${task.cpus} -a XX -a2 XX --ultra
        
    if [ ! -f ultraplex_demux_5bc_no_match_Rev.fastq.gz ]
    then
        touch ultraplex_demux_5bc_no_match_Rev.fastq.gz
        touch ultraplex_demux_5bc_no_match_Fwd.fastq.gz
    fi
    """
}
process trim_barcode_r1 {
    publishDir "${params.outDir}", mode: 'copy', overwrite: false
    input:
        val libName
        tuple val(rbp_label), val(rbp)
    output:
        path("ultraplex_demux_${rbp_label}_{Fwd,Rev}.fastq.Tr.fastq.gz"), optional: true, emit: trimmed_reads
        path("${rbp_label}_trim_barcode_r1.metrics")
    shell:
    """
    # rbptrimmed=\$(echo -n ${rbp_label} | tail -c +2 | head -c -1)
    bcstring=\$(grep -P ${rbp_label}\$ ${params.barcode_csv} | cut -f1 -d:)
    if [ -z "\$bcstring" ]
    then
        echo ${rbp[0]} ${rbp[1]} \$bcstring "Nothing Nothing"> ${rbp_label}_trim_barcode_r1.metrics
    else
        rev_bar=\$(echo \$bcstring | tr ACGTacgt TGCAtgca | rev)
        echo ${rbp[0]} ${rbp[1]} \$bcstring \$rev_bar > ${rbp_label}_trim_barcode_r1.metrics
        fastp -i ${rbp[0]} -I ${rbp[1]} -o ${rbp[0].baseName}.Tr.fastq.gz -O ${rbp[1].baseName}.Tr.fastq.gz --adapter_sequence \$rev_bar 2>> ${rbp_label}_trim_barcode_r1.metrics
        zcat ${rbp[1]} | grep -v "@" | grep \$rev_bar | wc -l >> ${rbp_label}_trim_barcode_r1.metrics
        zcat ${rbp[1].baseName}.Tr.fastq.gz | grep -v "@" | grep \$rev_bar | wc -l >> ${rbp_label}_trim_barcode_r1.metrics
    fi
    """
}
process fastqc_post_trim {
    publishDir "${params.outDir}", mode: 'copy', overwrite: false
    input:
        tuple val(rbp_label), path(rbp_path)
    output:
        path("ultraplex_demux_${rbp_label}_{Fwd,Rev}.fastq.Tr_fastqc.html")
        path("ultraplex_demux_${rbp_label}_{Fwd,Rev}.fastq.Tr_fastqc.zip")
    shell:
    """
    fastqc ultraplex_demux_${rbp_label}_Fwd.fastq.Tr.fastq.gz -t ${task.cpus}
    fastqc ultraplex_demux_${rbp_label}_Rev.fastq.Tr.fastq.gz -t ${task.cpus}
    """
}
process align_reads {
    cpus 8
    publishDir "${params.outDir}", mode: 'copy', overwrite: false
    input:
        tuple val(rbp_label), path(rbp_path)
    output:
        path("${rbp_label}.Aligned.sortedByCoord.out.bam"), emit: bam
    shell:
    """
    STAR \
        --alignEndsType EndToEnd \
        --genomeDir ${params.star_sjdb} \
        --genomeLoad NoSharedMemory \
        --outBAMcompression 10 \
        --outFileNamePrefix ${rbp_label}. \
        --winAnchorMultimapNmax 100 \
        --outFilterMultimapNmax 100 \
        --outFilterMultimapScoreRange 1 \
        --outSAMmultNmax 1 \
        --outMultimapperOrder Random \
        --outFilterScoreMin 10 \
        --outFilterType BySJout \
        --limitOutSJcollapsed 5000000 \
        --outReadsUnmapped Fastx \
        --outSAMattrRGline ID:${rbp_label} \
        --outSAMattributes All \
        --outSAMmode Full \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outStd Log \
        --readFilesIn ${rbp_path[0]} ${rbp_path[1]} \
        --readFilesCommand zcat \
        --runMode alignReads \
        --runThreadN ${task.cpus}
    """
}
process umi_dedup {
    tag "$rbp_label"
    cpus 8
    publishDir "${params.outDir}", mode: 'copy', overwrite: false
    input:
        path(bam)
    output:
        path("*.Aligned.sortedByCoord.out.bam.bai")
        path("*.Aligned.sortedByCoord.out.rmDup.bam")
        path("*.Aligned.sortedByCoord.out.rmDup.bam.bai")
    shell:
    """
    samtools index ${bam};
    umicollapse bam -i ${bam} -o ${bam.baseName}.rmDup.bam --umi-sep : --two-pass --paired;
    samtools index ${bam.baseName}.rmDup.bam;
    """
}
workflow { 
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
        
    tile_adapter(
        params.adapter_fwd, 
        params.adapter_rev
    )
    trim_adapter(
        params.libName, 
        tile_adapter.out.adapter_fwd_txt, 
        tile_adapter.out.adapter_rev_txt, 
        read_pairs_ch
    )
    extract_umi_and_trim_polyG(
        params.libName, 
        trim_adapter.out.read1,
        trim_adapter.out.read2
    )
    trim_umi_from_read2(
        params.libName,
        extract_umi_and_trim_polyG.out.read1,
        extract_umi_and_trim_polyG.out.read2
    )
    demultiplex(
        params.libName,
        trim_umi_from_read2.out.read1,
        trim_umi_from_read2.out.read2
    )
        | flatten
        | map {tuple(it.baseName.replace("ultraplex_demux_", "").replace("_Fwd.fastq", "").replace("_Rev.fastq", ""), it) }
        | groupTuple
        | set { rbps_group }
    trim_barcode_r1(
        params.libName,
        rbps_group
    ).trimmed_reads
        | flatten
        | map {tuple(it.baseName.replace("ultraplex_demux_", "").replace("_Fwd.fastq.Tr.fastq", "").replace("_Rev.fastq.Tr.fastq", ""), it)}
        | groupTuple
        | set { trimmed_rbps_group }
    fastqc_post_trim(
        trimmed_rbps_group
    )
    align_reads(
        trimmed_rbps_group
    )
    umi_dedup(
        align_reads.out.bam
    )
} 

