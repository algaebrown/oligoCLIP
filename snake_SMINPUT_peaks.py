import pandas as pd
#snakemake -j 40 -s snake_SMINPUT_peaks.py --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn={params.cores} -q home-yeo" --configfile config/call_input_peaks/config.yaml --use-conda -n 


# fine INPUT bams
encode_df = pd.read_csv('/home/hsher/projects/oligo_results/ENCODE_reprocess.csv', index_col = 0)
encode_df['encode_id'] = encode_df['uid'].str.split('_', expand = True)[0]
unique_encode_input = encode_df.sort_values(by = 'uid').drop_duplicates(subset = ['encode_id']) 

nature2016_df = pd.read_csv('/home/hsher/projects/oligo_results/nature2016_hek_full.csv', index_col = 0)
unique_nature_input = nature2016_df.drop_duplicates(subset = ['RBP'])

all_inputs = pd.concat([unique_encode_input, unique_nature_input], axis = 0)

print('uids', all_inputs['uid'].unique())

SPECIES=config['SPECIES']
SCRIPT_PATH=config['SCRIPT_PATH']
WORKDIR=config['WORKDIR']
workdir: WORKDIR

module peak_call:
    snakefile:
        "rules/Snake_CLIPper.py"
    config: config

module peak_anno:
    snakefile:
        "rules/Snake_peakanno.py"
    config: config

rule all:
    input:
        expand("output/CLIPper/{sample_label}.peaks.filtered.bed", sample_label = all_inputs['uid'].tolist())
    output:
        "done"
    params:
        error_out_file = "error_files/all",
        run_time = 1,
        cores = "1",
        memory = "20",
        job_name = "all"
    shell:
        "echo $(date)  > {output};"

use rule clipper from peak_call as peak_call_input with:
    input:
        bam = lambda wildcards: all_inputs.loc[all_inputs['uid']==wildcards.sample_label, 'bam_control_0'].iloc[0],
        bai = lambda wildcards: all_inputs.loc[all_inputs['uid']==wildcards.sample_label, 'bam_control_0'].iloc[0]+'.bai'
    output:
        peak = "output/CLIPper/{sample_label}.peaks.bed"

use rule filter_clipper from peak_anno as filter_input_peaks with:
    input:
        "output/CLIPper/{sample_label}.peaks.bed"
    output:
        "output/CLIPper/{sample_label}.peaks.filtered.bed"