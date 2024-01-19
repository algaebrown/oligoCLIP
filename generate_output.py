def preprocess_outputs():
    ''' return preprocessing outputs'''
    outputs = expand("{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam.bai", libname = libnames, sample_label = rbps
    )+expand("{libname}/bw/COV/{sample_label}.{strand}.bw", libname = libnames, sample_label = rbps, strand = ['pos', 'neg']
    )+['QC/fastQC_basic_summary.csv',
        'QC/fastQC_passfail.csv',
        'QC/cutadapt_stat.csv',
        "QC/mapping_stats.csv",
        "QC/dup_level.csv",
        'QC/demux_read_count.txt',
        "QC/summary.csv"
        ]+expand("counts/genome/vectors/{libname}.{sample_label}.counts",
        libname = libnames, sample_label = rbps
    )+expand("QC/read_count/{libname}.{metric}.csv", libname = libnames, metric = ['region', 'genetype', 'cosine_similarity']
    )+expand("counts/repeats/tables/{repeat_type}/{experiment}.{sample_label}.tsv.gz", 
    experiment = experiments, sample_label = rbps, repeat_type = ['name', 'class', 'family']
    )+expand("counts/repeats/megatables/{repeat_type}/{libname}.tsv.gz", libname = libnames, repeat_type = ['name', 'class', 'family'])
    return outputs

def debugging_output():
    ''' blast unmapped reads/reads without barcode'''
    outputs = []
    if config['debug']:
        outputs += expand("QC/unmapped_blast_output/{libname}.{sample_label}.short.blast.tsv", 
        libname = libnames, sample_label = rbps)
        # TODO: these rules only work for PE pipeline
        # outputs += expand("QC/unmapped_blast_output/{libname}.{sample_label}.1.blast.tsv",
        # libname = libnames, sample_label = rbps)
        # outputs += expand("QC/nobarcode_blast_output/{libname}.blast.tsv",libname = libnames, sample_label = rbps)
    return outputs

def skipper_outputs():
    ''' generate skipper outputs, enriched windows, finemapped windows and motifs'''
    outputs = []
    if not singleplex:
        # normalize to internal libraries
        outputs+=expand("skipper/{bg_sample_label}/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
        libname = libnames,
        clip_sample_label = list(set(rbps)-set(config['AS_INPUT'])), # cannot call on itself
        bg_sample_label = config['AS_INPUT'] if config['AS_INPUT'] else []
        )+expand("skipper/{bg_sample_label}/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
        libname = libnames,
        sample_label = list(set(rbps)-set(config['AS_INPUT'])), 
        signal_type = ['CITS', 'COV'], # cannot call on itself
        bg_sample_label = config['AS_INPUT'] if config['AS_INPUT'] else []
        )+expand("skipper/{bg_sample_label}/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
        libname = libnames,
        sample_label = config['RBP_TO_RUN_MOTIF'],
        signal_type = ['CITS', 'COV'],
        bg_sample_label = config['AS_INPUT'] if config['AS_INPUT'] else []
        )
        # normalize to complementary control
        outputs += expand("skipper_CC/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
        libname = libnames,
        clip_sample_label = list(set(rbps)-set(config['AS_INPUT'])), # cannot call on itself
        )+expand("skipper_CC/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
        libname = libnames,
        sample_label = list(set(rbps)-set(config['AS_INPUT'])), 
        signal_type = ['CITS', 'COV'] # cannot call on itself
        )+expand("skipper_CC/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
        libname = libnames,
        sample_label = config['RBP_TO_RUN_MOTIF'],
        signal_type = ['CITS', 'COV']
        )
    # normalize to external bams
    if external_normalization:
        # normalize to external library
        outputs+=expand("skipper_external/{external_label}/enriched_windows/{libname}.{clip_sample_label}.enriched_windows.tsv.gz",
        external_label = list(external_normalization.keys()),
        libname = libnames,
        clip_sample_label = list(set(rbps)-set(config['AS_INPUT']))
        )
    return outputs

def beta_binom_mixture_outputs():
    ''' generate output for beta-binomial mixture'''
    outputs = []
    if not singleplex:
        # complementary control
        outputs+=expand("beta-mixture_CC/{libname}.{sample_label}.enriched_windows.tsv",
        libname = libnames,
        sample_label = list(set(rbps)-set(config['AS_INPUT'])),
        signal_type = ['CITS', 'COV'],
        strand = ['pos', 'neg']
        )+expand("beta-mixture_CC/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.{strand}.bw",
        libname = libnames,
        sample_label = list(set(rbps)-set(config['AS_INPUT'])),
        signal_type = ['CITS', 'COV'],
        strand = ['pos', 'neg']
        )+expand("beta-mixture_CC/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
        libname = libnames,
        sample_label = config['RBP_TO_RUN_MOTIF'],
        signal_type = ['CITS', 'COV'])
        
        # internal control
        outputs+=expand("beta-mixture/{bg_sample_label}/{libname}.{clip_sample_label}.enriched_windows.tsv",
        libname = libnames,
        clip_sample_label = list(set(rbps)-set(config['AS_INPUT'])),
        bg_sample_label = config['AS_INPUT'] if config['AS_INPUT'] else [] #TODO: IgG secondary analysis
        )+expand("beta-mixture/{bg_sample_label}/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.{strand}.bw",
        libname = libnames,
        sample_label = list(set(rbps)-set(config['AS_INPUT'])),
        bg_sample_label = config['AS_INPUT'] if config['AS_INPUT'] else [],
        signal_type = ['CITS', 'COV'],
        strand = ['pos', 'neg']
        )+expand("beta-mixture/{bg_sample_label}/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
        libname = libnames,
        bg_sample_label = config['AS_INPUT'] if config['AS_INPUT'] else [],
        sample_label = config['RBP_TO_RUN_MOTIF'],
        signal_type = ['CITS', 'COV'])

    if external_normalization:
        outputs += expand("beta-mixture_external/{external_label}/{libname}.{clip_sample_label}.enriched_windows.tsv",
        libname = libnames, clip_sample_label = rbps,
        external_label = list(external_normalization.keys())
        )+expand("beta-mixture_external/{external_label}/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.{strand}.bw",
        libname = libnames,
        sample_label = list(set(rbps)-set(config['AS_INPUT'])),
        signal_type = ['CITS', 'COV'],
        strand = ['pos', 'neg'],
        external_label = list(external_normalization.keys())
        )+expand("beta-mixture_external/{external_label}/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html",
        libname = libnames,
        sample_label = config['RBP_TO_RUN_MOTIF'],
        external_label = list(external_normalization.keys()),
        signal_type = ['CITS', 'COV'])
    # complementart control bigwigs
    if len(set(rbps)-set(config['AS_INPUT']))>2 and len(rbps)>1:
        outputs+=expand("{libname}/bw_bg/COV/{sample_label}.{strand}.bw", libname = libnames, sample_label = rbps, strand = ['pos', 'neg']
        )
    return outputs

def DMN_outputs():
    ''' generate output from Dirichlet Multinomial mixture!'''
    outputs = []
    if not singleplex:
        outputs+=expand("DMM/{libname}.mixture_weight.tsv", libname = libnames
        )+expand("DMM/homer/finemapped_results/{signal_type}/{libname}.{sample_label}/homerResults.html", libname = libnames,
        sample_label = config['RBP_TO_RUN_MOTIF'],
        signal_type = ['CITS', 'COV']
        )+expand("DMM/finemapping/mapped_sites/{signal_type}/{libname}.{sample_label}.finemapped_windows.bed.gz",
        libname = libnames,
        sample_label = rbps,
        signal_type = ['CITS', 'COV']
        )+expand('mask/{libname}.genome_mask.csv',
        libname = libnames,
        )+expand('mask/{libname}.repeat_mask.csv',
        libname = libnames,
        )
    
    return outputs

def clipper_outputs():
    ''' generate output from CLIPper'''
    outputs = []
    if not singleplex:
        # internal background
        outputs+=expand("CLIPper.{bg}/{libname}.{sample_label}.peaks.normed.compressed.annotate.bed",
        bg = config['AS_INPUT'] if config['AS_INPUT'] else [],
        sample_label = list(set(rbps)-set(config['AS_INPUT'])),
        libname = libnames
        )
        # complementary control
        outputs+=expand("CLIPper_CC/{libname}.{sample_label}.peaks.normed.compressed.annotate.bed",
        sample_label = list(set(rbps)-set(config['AS_INPUT'])),
        libname = libnames
        )+expand("CLIPper_CC/{libname}.{sample_label}.peaks.normed.compressed.motif.svg",
        sample_label = config['RBP_TO_RUN_MOTIF'],
        libname = libnames
        )
    # external
    if external_normalization:
        outputs+= expand("CLIPper-{external_label}/{libname}.{sample_label}.peaks.normed.compressed.annotate.bed",
            sample_label = list(set(rbps)-set(config['AS_INPUT'])),
            libname = libnames,
            external_label = list(external_normalization.keys())
            )
    return outputs

def comparison_outputs():
    outputs = expand("comparison/piranha/CC/{libname}.{sample_label}.bed",
        libname = libnames,
        sample_label =list(set(rbps)-set(config['AS_INPUT'])),
    )
    # )+expand("comparison/pureclip/{libname}.{sample_label}.bind.bed",
    #     libname = libnames,
    #     sample_label = list(set(rbps)-set(config['AS_INPUT']))
    # ) # very slow to run
    # )+expand("comparison/omniCLIP/output/{libname}.{sample_label}.omniclip_done.txt",
    #     libname = libnames,
    #     sample_label = list(set(rbps)-set(config['AS_INPUT']))
    # ) # it dies all the time
    return outputs

def get_output(clipper, skipper, comparison):
    output = preprocess_outputs()
    
    
    output += DMN_outputs() + beta_binom_mixture_outputs()
    if skipper:
        output += skipper_outputs()
    if clipper:
        output += clipper_outputs()
    output += debugging_output()
    if comparison:
        output += comparison_outputs()
    return output
