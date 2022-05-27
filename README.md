# oligoCLIP
pipeline for oligoCLIP. Goes from fastq -> bam.

# How to run.
- you need to install snakemake and mamba.
- fill out the config file. see [example](https://github.com/algaebrown/oligoCLIP/blob/master/eclipse_multi.yaml). In this file, you will need to provide `fastq_menifest`, where it points to 1 or multiple libraries with the *SAME* barcode set. Look at this [example file](https://github.com/algaebrown/oligoCLIP/blob/master/multiplex1.csv) to see how it is structures. You will also need a `barcode` [csv](https://github.com/algaebrown/oligoCLIP/blob/master/barcodes.csv) telling the software what barcode are for which RBP. The RBP names must be unique. The rest of the parameters are `ADAPTOR_PATH`, and adaptor, `STAR_REP` star indicies build from repbase, used to remove the repetitive sequences. `STAR_DIR` is the index built from genome sequences. Lastly is `TOOL_PATH`, which should be the directory to this repository.
- Then you run. `snakemake -s snakeCLIP.py -j 12 --keep-going --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --directory /home/hsher/scratch/ABC_reprocess/ --configfile eclipse_multi.yaml -n` 
