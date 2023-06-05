# oligoCLIP: Antibody barcoded eCLIP(ABC) processing pipeline from fastq.gz to windows and motifs (Nextflow)
- [original ABC paper](https://www.nature.com/articles/s41592-022-01708-8): use `snakeABC_SE.py`
- Yeolab paired-end protocol: use ``

# Installation
- You need to have Nextflow

# How to run.
```
NXF_OPTS="-Dleveldb.mmap=false"
NXF_SINGULARITY_CACHEDIR=singularity/
nextflow -C nextflow.config run modules/preprocess/preprocess_pe.nf -profile singularity -resume
```

# Config
- Example: `` 
## Basic Inputs:
### `MANIFEST`: 
- Example: 
...
- columns:
...
### `barcode_csv`: specifying barcode sequencing per Antibody/RBP
- Example: `config/barcode_csv/iter5.csv`
- Notebook to generate this file (Yeolab internal user): `utils/generate barcode-iter5.ipynb`
- delimiter: `:`
- columns:
    - 1st column: barcode sequence
        - YeoLab internal protocol: read 2 starts with this sequence. Double check with `zcat READ2.fastq.gz | grep BARCODE_SEQ`. This sequence is reverse complement to the adapter sequence (see notebook for detail)
        - ABC: read starts with this sequence.
    - 2nd column: Antibody/RBP name, Should not contain space, special characters such as #,%,*.

### OTHER


## Dependencies:


## Preprocessing options:


## Annotations:


# Output Files
## Trimmed fastqs, bams, bigwigs:


## QC: Quality control files. 


## Peaks, Binding Sites


### Mudskipper: Mixture Modelling in `beta-mixture_*` and `DMM`.
...

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
- `skipper.smk` is entirely stolen from skipper
- `normalization_DMN` contains Mudskipper code, which does mixture model/generative clustering.
- `clipper.smk` runs CLIPper and the chi-square things. Stolen from ENCODE pipeline.
- `analysis.smk` and `finemap.smk`: runs finemapping, motif detection from Skipper and MudSkipper
- rules for bigwig generation is in `MAKE_TRACK` repository.
- `merge_bw.smk` sums up bigwigs to make complementary control.




