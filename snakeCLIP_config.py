# Configuration file 

import os
import sys
import glob

# Input data
umi_pattern="CCCCCNNNNN"
cell_number=10 # RBP number

FASTQ="/projects/ps-yeolab5/hsher/ABC.fastq.gz"
name="ABC_CLIP"
barcode="/home/hsher/projects/oligoCLIP_pipe/bc.tsv"

adaptor="InvRiL19"

# Resources
CHROM_SIZES = '/projects/ps-yeolab3/eboyle/resources/hg38.chrom.sizes'
BLACKLIST = None


STAR_DIR='/projects/ps-yeolab4/software/eclip/0.7.0/examples/inputs/star_2_7_gencode29_sjdb'
STAR_REP='/projects/ps-yeolab4/software/eclip/0.7.0/examples/inputs/star_2_7_homo_sapiens_repbase_fixed_v2'

ADAPTOR_PATH='/projects/ps-yeolab4/software/eclip/0.7.0/examples/inputs/'

PYTHON3_PATH='/home/hsher/miniconda3/bin/python' # TODO: containerize
