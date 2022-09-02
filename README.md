# oligoCLIP
pipeline for oligoCLIP contains 2 stages. `Stage1: SnakeCLIP.py`: Goes from fastq -> bam. `Stage2: SnakePeakMain`: takes bam files for IP and background, perform clipper peak calls, normalization and annotation/motif calling.

# How to run.
- you need to install snakemake and mamba.

# Stage 1:
- fill out the config file. see [example](https://github.com/algaebrown/oligoCLIP/blob/master/eclipse_multi.yaml). In this file, you will need to provide `fastq_menifest`, where it points to 1 or multiple libraries with the *SAME* barcode set. Look at this [example file](https://github.com/algaebrown/oligoCLIP/blob/master/multiplex1.csv) to see how it is structures. You will also need a `barcode` [csv](https://github.com/algaebrown/oligoCLIP/blob/master/config/barcode_csv/barcodes.csv) telling the software what barcode are for which RBP. The RBP names must be unique. The rest of the parameters are `ADAPTOR_PATH`, and adaptor, `STAR_REP` star indicies build from repbase, used to remove the repetitive sequences. `STAR_DIR` is the index built from genome sequences. Lastly is `TOOL_PATH`, which should be the directory to this repository.
- Then you run. `snakemake -s snakeCLIP.py -j 12 --keep-going --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --directory /home/hsher/scratch/ABC_reprocess/ --configfile eclipse_multi.yaml -n` 

# Stage 2:
- fill out config file [example](https://github.com/algaebrown/oligoCLIP/blob/master/config/peak_call_config/snake_CLIPper_downsample.yaml)
- `SCRIPT_PATH`: path to the `scripts/` folder in oligoCLIP
- `SPECIES`: CLIPper `--species` flag. 
- `MANIFEST`: bam manifest specifying which IP to pair with which background. [example](https://github.com/algaebrown/oligoCLIP/blob/master/config/peak_call_csv/downsample_peakcall.csv). column `bam_0` contains the IP bam. `bam_control_0` contains the background bam. Can be RNA-seq, IgG library and such. `uid` is the unique ID for each IP-background pair. All bams must be indexed (ideally outputs from Stage 1)
- `WORKDIR`: working directory

- GTF: *.gtf.db annotation used by [annotator](https://github.com/byee4/annotator). 
- ANNOTATOR_SPECIES: species name used by [annotator](https://github.com/byee4/annotator).
- GENOMEFA: fasta file containing genome sequence. 
- CLIPper_pvalthes: pvalue threshold to filter CLIPper, used by motif calling. All peaks above this threshold will be used to find motif.
