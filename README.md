# oligoCLIP
pipeline for oligoCLIP contains 2 stages. `Stage1: SnakeCLIP.py`: Goes from fastq -> bam. `Stage2: SnakePeakMain`: takes bam files for IP and background, perform clipper peak calls, normalization and annotation/motif calling.

# How to run.
- you need to install snakemake and mamba.

# Stage 1:
- fill out the config file. see [example](https://github.com/algaebrown/oligoCLIP/blob/master/eclipse_multi.yaml). In this file, you will need to provide `fastq_menifest`, where it points to 1 or multiple libraries with the *SAME* barcode set. Look at this [example file](https://github.com/algaebrown/oligoCLIP/blob/master/multiplex1.csv) to see how it is structures. You will also need a `barcode` [csv](https://github.com/algaebrown/oligoCLIP/blob/master/barcodes.csv) telling the software what barcode are for which RBP. The RBP names must be unique. The rest of the parameters are `ADAPTOR_PATH`, and adaptor, `STAR_REP` star indicies build from repbase, used to remove the repetitive sequences. `STAR_DIR` is the index built from genome sequences. Lastly is `TOOL_PATH`, which should be the directory to this repository.
- Then you run. `snakemake -s snakeCLIP.py -j 12 --keep-going --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --directory /home/hsher/scratch/ABC_reprocess/ --configfile eclipse_multi.yaml -n` 

# Stage 2:
- fill out config file [example](https://github.com/algaebrown/oligoCLIP/blob/master/config/peak_call_config/snake_CLIPper_downsample.yaml)
- `SCRIPT_PATH`: path to the `scripts/` folder in oligoCLIP
- `SPECIES`: CLIPper `--species` flag. 
- `MANIFEST`: bam manifest specifying which IP to pair with which background. [example](https://github.com/algaebrown/oligoCLIP/blob/master/config/peak_call_csv/downsample_peakcall.csv). column `bam_0` contains the IP bam. `bam_control_0` contains the background bam. Can be RNA-seq, IgG library and such. `uid` is the unique ID for each IP-background pair. All bams must be indexed (ideally outputs from Stage 1)
- `WORKDIR`: working directory

- GTF: *.gtf.db annotation used by [annotator](https://github.com/YeoLab/annotator). 
- ANNOTATOR_SPECIES: species name used by [annotator](https://github.com/YeoLab/annotator).
- GENOMEFA: fasta file containing genome sequence. 
- CLIPper_pvalthes: pvalue threshold to filter CLIPper, used by motif calling. All peaks above this threshold will be used to find motif.

# Notes:
- To build the DAG: ```snakemake --snakefile snakeCLIP.py --dag --configfile eclipse_multi.yaml | dot -Tpdf > dag.pdf```
- run command example: 

```bash
snakemake \
--use-singularity \
--singularity-args="--bind /home/${USER} --bind /oasis --bind /projects" \
-s snakeCLIP.py \
-j 1 \
--configfile ./tests/ABC_multiplex/eclipse_multi.yaml \
--directory ./tests/ABC_multiplex/outputs/ \
--rerun-incomplete \
-p
```

- singularity images will by default be found in:

```bash
./tests/ABC_multiplex/outputs/.snakemake/singularity
```

But may be modified with:

```bash
--singularity-prefix
```

# Tests

- snakemake command must be run *without* ```--generate-unit-tests``` before running ```--generate-unit-tests``` (see above command)
- Once the original command is run, then you may generate the tests:

```bash
snakemake \
--use-singularity \
--singularity-args="--bind /home/${USER} --bind /oasis --bind /projects" \
-s snakeCLIP.py \
-j 1 \
--configfile ./tests/ABC_multiplex/eclipse_multi.yaml \
--directory ./tests/ABC_multiplex/outputs/ \
--generate-unit-tests \
--rerun-incomplete \
-p
```

- unit tests will by default be found in: 

```bash
./tests/ABC_multiplex/outputs/.tests/
```

- Few issues to be aware of:
  - generated pytests cannot recognize reading in "fastq_menifest" from config file. Possible solution is to include:
    - fastq_menifest contents inside config (eg.
    
  	  ```
        fastqs: 
          - libname: K562_rep1
            fastq: /home/bay001/projects/codebase/oligoCLIP/tests/ABC_multiplex/inputs/fastqs/ABC009_1_S1_R1_001_small.fastq.gz
          - libname: K562_rep2
            fastq: /home/bay001/projects/codebase/oligoCLIP/tests/ABC_multiplex/inputs/fastqs/ABC009_2_S2_R1_001_small.fastq.gz
          - libname: HEK293_rep1
            fastq: /home/bay001/projects/codebase/oligoCLIP/tests/ABC_multiplex/inputs/fastqs/ABC009_3_S3_R1_001_small.fastq.gz
          - libname: HEK293_rep2
            fastq: /home/bay001/projects/codebase/oligoCLIP/tests/ABC_multiplex/inputs/fastqs/ABC009_4_S4_R1_001_small.fastq.gz
  	  ```
      
    )
  - Singularity bind args not carried over into subprocess.check_output methods. Possible solutions include:
    - Ensuring singularity binds to the appropriate directories after each generated sp.check_output() command.
    - Run a sed replace command (eg. 
    
      ```
      sed -i 's/])/], env\=common.get_env())/g' *.py
      ```
    ). 
    **Note** This will also modify common.py too, which is distinct from test_* and should be reverted back to the original file (without calling get_env())
    
      ```common.get_env()```
    
      ```python
      def get_env():
          env = dict(os.environ, SINGULARITY_BINDPATH='/home,/projects,/oasis')
          return env
      ```
    

  - generate-unit-tests **do not properly escape config files**. [See issue](https://github.com/snakemake/snakemake/issues/843). Possible solutions include:
    - Manually going into each test and ensuring strings are quote-escaped (eg. 
        ```python
        ...
        "--configfile",
        "/projects/ps-yeolab3/bay001/codebase/oligoCLIP/tests/ABC_multiplex/eclipse_multi.yaml",
        ...
        ```
    )
  
  - generate-unit-tests **do not point to Snakefiles**
    - Manually ensure -s flag is set (eg. 
        ```python
        ...
        "-s",
        "/projects/ps-yeolab3/bay001/codebase/oligoCLIP/snakeCLIP.py",
        ...
        ```
    )

  - generate-unit-tests **do not properly format multi-file outputs**
    - Manually listing each output file as separate targets (eg.
    
        ```python
        "python",
        "-m",
        "snakemake", 
        "K562_rep1/bams/dros/IGF2BP2.Aligned.out.bam K562_rep1/bams/dros/IGF2BP2.Unmapped.out.mate1 K562_rep1/bams/dros/IGF2BP2.Log.final.out",
        ...
        ```
        <br>
        Should be:
        <br>

        ```python
        "python",
        "-m",
        "snakemake", 
        "K562_rep1/bams/dros/IGF2BP2.Aligned.out.bam",
        "K562_rep1/bams/dros/IGF2BP2.Unmapped.out.mate1",
        "K562_rep1/bams/dros/IGF2BP2.Log.final.out",
        ...
        ```
    
    )
  - generate-unit-tests **will only perform byte-by-byte comparisons**
    - This is problematic whenever files contain non-reproducible outputs (eg. timestamps). To fix this, within each unit test, you will need to extend ```OutputChecker()``` inside ```common.py``` so that each output can be properly compared. See ```class BamOutputChecker``` which overwites byte-checking BAM files (as BAM file headers differ) to instead check that the number of mapped reads is identical.

- Run a pytest to ensure all tests pass. Inside ```outputs/```, use this command: 

```bash
python -m pytest .tests/unit/ > log.txt 2>&1
```
  
- log.txt contents:

```
============================= test session starts ==============================
platform linux -- Python 3.10.5, pytest-7.1.2, pluggy-1.0.0
rootdir: /projects/ps-yeolab3/bay001/codebase/oligoCLIP/tests/ABC_multiplex/outputs
collected 21 items

.tests/unit/test_align_reads_to_Drosophila.py .                          [  4%]
.tests/unit/test_align_reads_to_REPEAT.py .                              [  9%]
.tests/unit/test_align_to_GENOME.py .                                    [ 14%]
.tests/unit/test_all.py .                                                [ 19%]
.tests/unit/test_cutadapt_round_one.py .                                 [ 23%]
.tests/unit/test_cutadapt_round_two.py .                                 [ 28%]
.tests/unit/test_demultiplex.py .                                        [ 33%]
.tests/unit/test_extract_umi.py .                                        [ 38%]
.tests/unit/test_fastqc_post_trim.py .                                   [ 42%]
.tests/unit/test_index_genome_bams.py .                                  [ 47%]
.tests/unit/test_make_good_barcode_tsv.py .                              [ 52%]
.tests/unit/test_mapstat_gather_dros.py .                                [ 57%]
.tests/unit/test_mapstat_gather_genome.py .                              [ 61%]
.tests/unit/test_mapstat_gather_repeat.py .                              [ 66%]
.tests/unit/test_qc_fastqc.py .                                          [ 71%]
.tests/unit/test_qc_trim1.py .                                           [ 76%]
.tests/unit/test_qc_trim2.py .                                           [ 80%]
.tests/unit/test_remove_barcode_and_reverse_complement.py .              [ 85%]
.tests/unit/test_sort_and_gzip_fastq.py .                                [ 90%]
.tests/unit/test_sort_bams.py .                                          [ 95%]
.tests/unit/test_umi_dedup.py .                                          [100%]

======================== 21 passed in 125.27s (0:02:05) ========================
```
