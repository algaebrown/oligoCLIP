#!/usr/bin/env nextflow

params.tile_length = 10
params.adapter_fwd = "AGATCGGAAGAGCACACGTC"
params.adapter_rev = "AGATCGGAAGAGCGTCGTGT"
params.main_path = "/projects/ps-yeolab3/bay001/codebase/oligoCLIP/"
params.quality_cutoff = 6
params.libName = 'testing'
params.umi_length = 10

params.greeting = 'Preprocessing eCLIP data!' 

params.reads = "/home/bay001/projects/encode5/temporary_data/merged_fastqs/small_test_R{1,2}.fastq.gz"
reads_ch = Channel.fromFilePairs(params.reads) 

greeting_ch = Channel.of(params.greeting) 

process tile_adapter {
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
    input:
        val libName
        val adapter_fwd_txt
        val adapter_rev_txt
        tuple val(replicateId), path(reads)
    output:
        path "${libName}/all.Tr.R1.fastq.gz"
        path "${libName}/all.Tr.R2.fastq.gz"
        /* path "${libName}.Tr.metrics" */
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
        -o ${libName}/all.Tr.R1.fastq.gz \
        -p ${libName}/all.Tr.R2.fastq.gz \
        --cores=0 \
        ${reads[0]} ${reads[1]} > ${libName}.Tr.metrics
        """
}
process extract_umi_and_trim_polyG {
    input:
        val libName
        path "${libName}/all.Tr.R1.fastq.gz"
        path "${libName}/all.Tr.R2.fastq.gz"
    output:
        path "${libName}/all.Tr.umi.R1.fastq.gz"
        path "${libName}/all.Tr.umi.R2.fastq.gz"
        path "QC/${libName}.umi.json"
        path "QC/${libName}.umi.html"

    script:
    """
    mkdir QC;
    fastp -i ${libName}/all.Tr.R1.fastq.gz \
        -I ${libName}/all.Tr.R2.fastq.gz \
        -o ${libName}/all.Tr.umi.R1.fastq.gz \
        -O ${libName}/all.Tr.umi.R2.fastq.gz \
        --disable_adapter_trimming \
        --umi \
        --umi_len=${params.umi_length} \
        --umi_loc=read1 \
        --trim_poly_g \
        -j QC/${libName}.umi.json \
        -h QC/${libName}.umi.html \
        -w ${task.cpus}
    """
}

workflow { 
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    tile_adapter(params.adapter_fwd, params.adapter_rev)
    trim_ch = trim_adapter(
        params.libName, 
        tile_adapter.out.adapter_fwd_txt, 
        tile_adapter.out.adapter_rev_txt, 
        read_pairs_ch
    )
    extract_umi_ch = extract_umi_and_trim_polyG(
        params.libName, 
        trim_ch
    )
    
} 