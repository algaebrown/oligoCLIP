#!/usr/bin/env nextflow

params.tile_length = 10
params.adapter_fwd = "AGATCGGAAGAGCACACGTC"
params.adapter_rev = "AGATCGGAAGAGCGTCGTGT"
params.main_path = "/projects/ps-yeolab3/bay001/codebase/oligoCLIP/"

params.greeting = 'Preprocessing eCLIP data!' 
greeting_ch = Channel.of(params.greeting) 

process tile_adaptor {
    input:
        val adapter_fwd
        val adapter_rev
    output:
        path "fwd_adapters.txt"
        path "rev_adapters.txt"
    script:
        """
        python ${params.main_path}/scripts/create_adaptor_tile.py ${params.adapter_fwd} fwd_adapters.txt ${params.tile_length}
        python ${params.main_path}/scripts/create_adaptor_tile.py ${params.adapter_rev} rev_adapters.txt ${params.tile_length}
        """
}

workflow { 
    tiles = tile_adaptor(params.adapter_fwd, params.adapter_rev)
} 