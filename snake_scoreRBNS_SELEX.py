# This file calculate motif from from RBNS motif for RBFOX2.
RBPNAME='RBFOX2'
rule get_sequence:
    input:
        "{something}.bed"
    output:
        "{something}.fa"
    params:
        run_time = 1,
        cores = "1",
        fa=config['GENOMEFA']
    shell:
        '''
        module load bedtools
        bedtools getfasta -fo {output} -fi {params.fa} -bed {input} -s
        '''
rule get_RBNS_SELEX_score:
    input:
        "{something}.fa"
    output:
        "{something}.motifscore.csv"
    params:
        run_time = "2",
        cores = 1,
        RBPname=RBPNAME,
        script_path=config['SCRIPT_PATH']
    conda:
        "envs/metadensity.yaml"
    shell:
        '''
        python {params.script_path}score_by_RBNS_SELEX.py {params.RBPname} {input} {output}
        '''