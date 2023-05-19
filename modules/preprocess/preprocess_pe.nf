#!/usr/bin/env nextflow

params.tile_length = 10
params.adapter_fwd = "AGATCGGAAGAGCACACGTC"
params.adapter_rev = "AGATCGGAAGAGCGTCGTGT"
params.main_path = "/projects/ps-yeolab3/bay001/codebase/oligoCLIP/"
params.quality_cutoff = 6
params.libname = 'testing'
params.fq1 = '/home/bay001/projects/encode5/temporary_data/merged_fastqs/small_test_R1.fastq.gz'
params.fq2 = '/home/bay001/projects/encode5/temporary_data/merged_fastqs/small_test_R2.fastq.gz'

params.greeting = 'Preprocessing eCLIP data!' 
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
        val libname
        val fq1
        val fq2
        val adapter_fwd_txt
        val adapter_rev_txt
    output:
        path "${libname}/all.Tr.fq1.gz"
        path "${libname}/all.Tr.fq2.gz"
        path "${libname}.Tr.metrics"
    script:
        """
        mkdir ${libname};
        cutadapt -a file:${adapter_fwd_txt} \
        -A file:${adapter_rev_txt} \
        --times 2 \
        -e 0.1 \
        --quality-cutoff ${params.quality_cutoff} \
        -m 23 \
        -o ${libname}/all.Tr.fq1.gz \
        -p ${libname}/all.Tr.fq2.gz \
        --cores=0 \
        ${fq1} ${fq2} > ${libname}.Tr.metrics
        """
}
workflow { 
    tile_adapter(params.adapter_fwd, params.adapter_rev)
    trimmed = trim_adapter(params.libname, params.fq1, params.fq2, tile_adapter.out.adapter_fwd_txt, tile_adapter.out.adapter_rev_txt)
} 