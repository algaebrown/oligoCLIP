

# This file does peak call, normalization, annotation and calculate RBFOX2 motif score
#snakemake -j 16 -s Snake_PeakMain.py --use-conda --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn={params.cores} -q home-yeo" --configfile config/snake_CLIPper_downsample.yaml --keep-going output/CLIPper/{Dan_singleplex_HEK293_rep1_RBFOX2,Dan_singleplex_HEK293_rep2_RBFOX2,ENCODE_Dan_singleplex_HEK293_rep1_RBFOX2,ENCODE_Dan_singleplex_HEK293_rep2_RBFOX2}.peaks.motifscore.csv
#snakemake -j 12 -s Snake_PeakMain.py --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn={params.cores} -q home-yeo" --configfile config/snake_CLIPper_encode.yaml -n 
#snakemake -j 12 -s Snake_PeakMain.py --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn={params.cores} -q home-yeo" --configfile config/snake_CLIPper_nature2016.yaml -n 
#snakemake -j 12 -s Snake_PeakMain.py --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn={params.cores} -q home-yeo" --configfile config/snake_CLIPper_katie.yaml --keep-going output/CLIPper/{katieoligo_RBFOX2_rep1,katieoligo_RBFOX2_rep2,katieoligo_RBFOX2_rep3,katieoligo_RBFOX2_rep4}.peaks.motifscore.csv output/{katieoligo_RBFOX2_rep1,katieoligo_RBFOX2_rep2,katieoligo_RBFOX2_rep3,katieoligo_RBFOX2_rep4}.peaks.normed.compressed.motifscore.csv
#snakemake -j 40 -s Snake_PeakMain.py --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn={params.cores} -q home-yeo" --configfile config/peak_call_config/snake_CLIPper_downsample.chi.yaml --use-conda -n 

from importlib.resources import path
import pandas as pd
import os



try:
    MANIFEST=config['MANIFEST']
    SPECIES=config['SPECIES']
    SCRIPT_PATH=config['SCRIPT_PATH']
    WORKDIR=config['WORKDIR']
    workdir: WORKDIR

    manifest = pd.read_table(MANIFEST, index_col = False, sep = ',')
    sample_labels = manifest.uid.tolist()
    all_rbfox = [s for s in sample_labels if 'RBFOX2' in s or '676' in s]
except Exception as e:
    print('config:', e)

module peak_anno:
    snakefile:
        "rules/annotate_peak.py"
    config: config

module peak_call:
    snakefile:
        "rules/clipper_and_norm.py"
    config: config

module motif:
    snakefile:
        "rules/snake_scoreRBNS_SELEX.py"
    config: config

module chi:
    snakefile:
        "rules/complementary_control.py"
    config: config

rule all:
    input:
        expand("output/CLIPper/{sample_label}.peaks.bed", sample_label = sample_labels)+
        expand("output/{sample_label}.peaks.normed.compressed.bed", sample_label = sample_labels),
        expand("output/{sample_label}.peaks.normed.compressed.annotate.bed", sample_label = sample_labels),
        expand("output/motif/{sample_label}.peaks.normed.compressed.filtered.annotate.svg", sample_label = sample_labels),
        expand("output/CLIPper/{sample_label}.peaks.annotate.bed", sample_label = sample_labels),
        expand("output/CLIPper/{sample_label}.peaks.filtered.svg", sample_label = sample_labels),
        expand("output/{sample_label}.peaks.normed.compressed.motifscore.csv", sample_label = all_rbfox),
        expand("output/CLIPper/{sample_label}.peaks.motifscore.csv", sample_label = all_rbfox),
        expand("summary/normalized/{sample_label}.peaks.summary", sample_label = sample_labels),
        expand("summary/CLIPper/{sample_label}.peaks.summary", sample_label = sample_labels),
        expand("chisq/motif/{combinations}.peaks.normed.compressed.filtered.annotate.svg", combinations = (manifest['lib']+'/'+manifest['uid']).tolist()),
        expand("chisq/{combinations}.peaks.summary", combinations = (manifest['lib']+'/'+manifest['uid']).tolist())

    output:
        "snakeCLIP.txt"
    params:
        error_out_file = "error_files/all",
        run_time = 1,
        cores = "1",
        memory = "20",
        job_name = "all"
    shell:
        """
        echo $(date)  > {output};
        echo created by HLH and the Yeo lab >> {output}
        """

use rule * from peak_call as call_*

use rule * from peak_anno as anno_*

use rule * from motif as motif_*

use rule * from chi as chi_*