# oligoCLIP: Antibody barcoded eCLIP(ABC) processing pipeline from fastq.gz to windows and motifs
- [original ABC paper](https://www.nature.com/articles/s41592-022-01708-8): use `snakeABC_SE.py`
- Yeolab paired-end protocol: use `snakeOligoCLIP_PE.py`

# Installation
- You need to have Snakemake:
    - snakemake instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
    - Yeolab internal users: `module load snakemake`
- Download this repository: `https://github.com/algaebrown/oligoCLIP.git` and then checkout this branch `git checkout oligo-pe`
    - `git branch` to double check you are on the right branch.
- Download depending repository and modify config variables as follow: # TODO: containerize or make to snakemake hub
    - Yeolab internal users don't need to.
- Install skipper dependecies and modify the following config variables:`JAVA_PATH`,`UMICOLLAPSE_PATH`, `R_EXE`. # TODO: containerize
    - follow [skipper instructions](https://github.com/YeoLab/skipper#prerequisites) to set up
- Most dependencies are already specified in `rules/envs`. When running snakemake, using `--use-conda` should automatically install everything for you.


# How to run.
1. prepare `PATH_TO_YOUR_CONFIG`. See below and `config/preprocess_config/oligope_iter5.yaml`
2. Run snakemake
```
snakemake -s snakeOligoCLIP_PE.py \ 
    -j 12 \
    --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" \
    --configfile PATH_TO_YOUR_CONFIG \
    --use-conda \
    --conda-prefix /home/hsher/snakeconda -npk
```
- `-s`: use `snakeOligoCLIP_PE.py` if you did YeoLab internal pair-end protocol. use `snakeABC_SE.py` if you did ABC
- `--configfile`: yaml file to specify your inputs, including where are the fastqs, what are the barcode, what reference genome...etc.
- the rest just snakemake command line options. [see documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
- `-j`: number of jobs to run at a same time
- `--cluster`: command to submit jobs to cluster. 
- `--use-conda`: ask snakemake to install everything for you using conda
- `--conda-prefix`: specify a fixed location to store conda envs to prevent snakemake installing them multiple times
- `-n`: dry run.
- `-k`: keep going even if something failed
- `-p`: print out command

# Config
- Example: `config/preprocess_config/oligope_iter5.yaml` 
## Basic Inputs:
### `MANIFEST`: a csv specifying fastq locations, replicates
- Example: 
    - YeoLab internal: `config/fastq_csv/katie_pe_iteration5.csv`
    - ABC: `config/fastq_csv/ABC_2rep.csv`
- columns:
    - `fastq1`&`fastq2`: *.fastq.gz file for read1 and read 2
    - `libname`: unique names for each library. Should not contain space, special characters such as #,%,*
    - `experiment`: unique names for experiment. **Rows with the same `experiment` will be treated as replicates.** Should not contain space, special characters such as #,%,*
### `barcode_csv`: specifying barcode sequencing per Antibody/RBP
- Example: `config/barcode_csv/iter5.csv`
- Notebook to generate this file (Yeolab internal user): `utils/generate barcode-iter5.ipynb`
- delimiter: `:`
- columns:
    - 1st column: barcode sequence
        - YeoLab internal protocol: read 2 starts with this sequence. Double check with `zcat READ2.fastq.gz | grep BARCODE_SEQ`. This sequence is reverse complement to the adapter sequence (see notebook for detail)
        - ABC: read starts with this sequence.
    - 2nd column: Antibody/RBP name, Should not contain space, special characters such as #,%,*.

### Outputs
- `WORKDIR`: output directory
- `RBP_TO_RUN_MOTIF`: list of RBP names to run motif analysis. Must be one of the rows in `barcode_csv`.
- `run_clipper`: True if you want CLIPper outputs (works, but slow)
- `run_skipper`: True if you want to run Skipper. (usually doesn't work in ABC)
- `run_comparison`: True if you want to run Piranha
- debug: True if you want to debug. This tries to blast the unmapped reads.

### Choosing backgrounds
By default if the below are left blank, we run Dirichlet Multinomial Mixture(DMM) for multiplex datasets, where RBPs are explicitly compared with each other. DMM is the best model for multiplex dataset. 

Unfortunately, DMM doesn't work for singleplex. Calling singleplex binding sites require "external control" (see below). Otherwise it will just stop at the read counting stage.

But if you want to add an background library, here is how to do:
#### "Internal control": a barcode that measures the background. They are in the same `fastq.gz`
- `AS_INPUT`: if you have a IgG antibody that everything will normalize against, type its name here. Must be one of the rows in `barcode_csv`. This can the background for skipper, CLIPper, and beta-binomial mixture model
#### "External control": a library that is NOT in the same fastq as your oligoCLIP/ABC
- specify them in `external_bam`
- This can be an eCLIP SMInput, total RNA-seq, IgG pull down from another experiment, bead control, spike-ins
- these will also be used as a background in skipper, CLIPper and beta-binomial mixture model
- the bams must be processed with the exact same STAR index as `STAR_DIR`, and is recommended to be processed with the same/similar mapping parameters as this repo or skipper.



## Dependencies:
- `SCRIPT_PATH`: Absolute path to `scripts` folder.
- `JAVA_PATH`,`UMICOLLAPSE_PATH`, `R_EXE`: skipper dependencies. See `Installation`.

## Preprocessing options:
- `adaptor_fwd`,`adaptor_rev`: adapter sequence to trim. Do not include barcode
- `tile_length`: we tile adapter sequences of this length so that indels don't mess up with trimming
- `QUALITY_CUTOFF`: default 15. cutadapt params
- `umi_length`: Length of unique molecular identifier (UMI).
- `STAR_DIR`: directory to STAR index

## Annotations:
- skipper annotations: [follow skipper instructions](https://github.com/YeoLab/skipper#prerequisites) or generate with [skipper_utils](https://github.com/algaebrown/skipper_utils)
    - Yeolab internal users: Brian had all sorts of annotations here `/projects/ps-yeolab4/software/skipper/1.0.0/bin/skipper/annotations/`.
- `CHROM_SIZES`
- `GENOMEFA`

# Output Files
## Trimmed fastqs, bams, bigwigs:
These are in the `EXPERIMENT_NAME` folders. For example, in your manifest.csv, there are two experiments, "GN_1019" and "GN_1020", then, under the `GN_1019/` folder you would see the following:
1. `fastqs`: The trimmed and the demultiplexed fastqs.
    - `all.Tr.fq*` is the adapter and UMI removed ones
    - `ultraplex*RBP.fastq.gz` are the demultiplexed.
2. `bams`: bam files.
    - `*rmDup.Aligned.sortedByCoord.out.bam` is the PCR deduplicated bams. This is usually the bam you would want to use for other analysis!
    - `*Aligned.sortedByCoord.out.bam` is before deduplication.
3. `bw` contains bigwig files!
    - `CITS` calculates the number of 5' read stops per nucleotide position. 5' read stop is presumed to be the crosslinking site.
    - `COV` calculates read coverage.
4. `bw_bg` contains bigwigs of "complementary control(CC)", namely, summing all the other RBPs together.

## QC: Quality control files. 
These contain numbers that will help you debug your protocol. To look at them swiftly, I personally like to install VS code, and install some csv editing extension so the file can appear like a table. Everything is compiled into `QC/summary.csv`.
1. `cutadapt_stat.csv` contains how many times reads contain adapter, and how many reads are too short after adapter trimming.
    - One of my favorite column is "Reads/Pairs that were too short". When reads are too short after trimming, it means they are probably adapter dimer, or you fragment it too much.
2. `fastQC*` are the summary from fastqs that has been trimmed. 
    - CLIP reads are usually bad at GC content. So it is normal to see failed GC.
    - I will always look at "adapter content" in `fastQC_passfail.csv` to make sure adapters are all gone. If this column fails, maybe you input the wrong adapter sequence, or gets contaminated by some other stuffs.
3. Read count after demultiplexing. `demux_read_count.txt`.
4. Duplication: `dup_level.csv`. This file contains the number of reads after and before PCR deduplication. If the column "unique reads" is very low, you might have amplified too much!
5. Mapping statistics `mapping_stat.csv`.
    - "Unique mapped reads%" are the percentage of reads that map to only 1 position in the genome. These were called "Usable reads" in the old eCLIP terminology.
    - "% of reads mapped to multiple loci" are the reads that map to 2~100 positions. These are mostly ribosomal RNAs. It is normal to see quite some (30-50%) multimapping reads in eCLIP
    - Then it is the reason why the rest are unmapped:
        - "% of reads mapped to multiple loci", "% of reads unmapped: too short", and "other reasons". These will help your bioinformatician debug what is going wrong. Common reasons to have lots of unmapped reads:
            - Adapter trimming is bad. The read still contains adapter sequences when they enter mapping algorithms.
            - Cells is contaminated. The baterial/fungal genomes gets sequences and does not map to human genome.
6. Read level metrices in `read_count/` folder:
    - `*cosine_similarity.csv`: Here we construct a vector, contains read counts per genomic bin for each RBP. Cosine similarity measure globally how two RBP are similar with each other. If you see big numbers, it means the two RBPs are very similar. In the Antibody barcoded CLIP paper we published, the cosine similarity is typically between 0.4~0.7. If all of the RBPs are similar, it can indicate loss of specificity.
    - `*gene_name.csv`: This contains how many read per gene for each RBP.
    - `*gene_type.csv`: contains how many read per gene/transcript type for each RBP.
    - `*region.csv`: This contains how many read per region for each RBP.
    - All of the above can help you see if there is the right enrichment and specificity for your multiplex CLIP.

## Peaks, Binding Sites
This pipeline tries to integrate multiple peak callers/binding site finders and orchestrate secondary analysis (peaks per region, motifs etc).

### Mudskipper: Mixture Modelling in `beta-mixture_*` and `DMM`.
- `beta-mixture*` considers enrichment of RBP reads against complementary control (CC) or an internal library(IgG)! It is not great at detecting shared binding sites. But it runs fast.
    - For biologist, all that matters is `*enriched_windows.tsv`. This contains the presumed binding sites.
        - column `logLR` measures confidence. The number represents the log likelihood ratio(LR) of a window being a binding site versus not, aka, how likely is it to observe the data if it is bound, versus it being not bound. A logLR of 2 means about 10x more likely to observe the data from binding than being something sh**ty. Higher means high confidence of it being a binding site
        - `p_bar`,`p_raw`, `fc_raw`, `fc_bar`: Effect sizes
            - `*_bar` is the estimate from the model.
            - `*_raw` is directly calculated from counts:
            - p= (read in RBP)/(read in window)
            - fc = ((read in RBP)/(read in window))/((all reads in RBP)/total reads)
    - Secondary analysis: 
        - `feature_logistic.pdf` and `feature_ridge.pdf` contains how likely each types of regions is bound.
        - `*summary.csv`: What gene types/transcript types/region types are bound.
        - `homer/` folder contains motifs
        - `finemapping` contains finemapped binding sites. It isn't great for splice site binding proteins. TODO: make algorithm better,
    - Modelling outputs:
        - `fit.rda` contains everything including models of various numbers of components(K).
        - `*alpha.tsv` contains the parameter alpha/beta for beta-binomial distribution for each component.
        - `*null_alpha.tsv` contains the parameter alpha/beta for a single component beta-binomial distribution.
        - `*goodness_of_fit.pdf`: AIC, BIC, log likelihood for model selection.
        - `*label_component.csv`: The labels (bound vs not bound) for each component.
        - `*mixture_weight.tsv`: This is namely badly. What is contains is $E[z_ik]$ which is "how likely each window belong to a cluster".
        - `*weights`: The "mixture weight".
- `DMM`(Dirichlet Multinomial Mixture) considers the distribution of reads among all RBPs without summing the rest into CC. This model detects shared binding site better than beta-mixture model. But it is so f* slow. (but still faster than CLIPper. ^.<)
    - The folder structure is the same.
    - Biologists would only care for `*enriched_windows.tsv` which are the binding sites.

### Skipper: in `skipper*`
- `skipper_CC` models RBP versus complementary control.
- `skipper_{INTERNAL_CONTROL_NAME}` models RBP versus an internal control library, e.g. IgG.
- `skipper_external/{EXTERNAL_CONTROL_NAME}` models RBP versus an external control library, e.g. RNA-seq or SMInput.
- For skipper outputs, see skipper's documentation!

### CLIPper: in `CLIPper*`
- `CLIPper` only uses the IP to find local read enrichment.
- `CLIPper_CC` contains local read enrichment, "normalized to" complementary control using chi-square or fisher exact test. This is what we publised in the original paper.
- `CLIPper-{EXTERNAL_CONTROL_NAME}`: contains peaks "normalized to" an external control library. (SMInput or total RNA-seq)
- `CLIPper.{INTERNAL_CONTROL_NAME}`: contains peaks "normalized to" an internal control library. (IgG)
- `*normed.compressed.annotate.bed` is the final output. See the ENCODE pipeline for columns specification

### Comparison: `comparison/` Only if you want to run Piranha, OmniCLIP and PureCLIP.
- See their respective documentation for detail.


# For developers:
The meat of the pipeline is in `rules/`:
- `se_preprocess` and `pe_preprocess` takes fastq --> trim -> demultiplex -> deduplicate -> bams
- `QC` contains rules to assemble quality control statistics, and some additional debugging rules such as investigating unmapped reads and those without barcode.
- rules for bigwig generation is in `make_track.smk`
- `merge_bw.smk` sums up bigwigs to make complementary control.
- `normalization_DMN`, `repeat_DMN` contains Mudskipper code, which does mixture model/generative clustering in genomic windows and repeat windows.
- `skipper.smk`, `repeat.smk` is entirely stolen from skipper
- `clipper.smk` runs CLIPper and the chi-square things. Stolen from ENCODE pipeline.
- `analysis.smk` and `finemap.smk`: runs finemapping, motif detection from Skipper and MudSkipper

Some rules to help you debug
- `map_r1.smk` and `multimap.smk`





