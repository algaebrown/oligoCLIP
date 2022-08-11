import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_demultiplex():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/demultiplex/data")
        expected_path = PurePosixPath(".tests/unit/demultiplex/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("K562_rep1/fastqs/IGF2BP2.umi.fqTrTr.fastq K562_rep1/fastqs/RBFOX2.umi.fqTrTr.fastq K562_rep1/fastqs/PUM2.umi.fqTrTr.fastq K562_rep1/fastqs/FAM120A.umi.fqTrTr.fastq K562_rep1/fastqs/DDX3.umi.fqTrTr.fastq K562_rep1/fastqs/ZC3H11A.umi.fqTrTr.fastq K562_rep1/fastqs/EIF3G.umi.fqTrTr.fastq K562_rep1/fastqs/PRPF8.umi.fqTrTr.fastq K562_rep1/fastqs/LIN28B.umi.fqTrTr.fastq K562_rep1/fastqs/SF3B4.umi.fqTrTr.fastq K562_rep1/barcode.log", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "K562_rep1/fastqs/IGF2BP2.umi.fqTrTr.fastq",
            "K562_rep1/fastqs/RBFOX2.umi.fqTrTr.fastq",
            "K562_rep1/fastqs/PUM2.umi.fqTrTr.fastq",
            "K562_rep1/fastqs/FAM120A.umi.fqTrTr.fastq",
            "K562_rep1/fastqs/DDX3.umi.fqTrTr.fastq",
            "K562_rep1/fastqs/ZC3H11A.umi.fqTrTr.fastq",
            "K562_rep1/fastqs/EIF3G.umi.fqTrTr.fastq",
            "K562_rep1/fastqs/PRPF8.umi.fqTrTr.fastq",
            "K562_rep1/fastqs/LIN28B.umi.fqTrTr.fastq",
            "K562_rep1/fastqs/SF3B4.umi.fqTrTr.fastq",
            "K562_rep1/barcode.log",
            "-f", 
            "-j1",
            "--keep-target-files",
            "--configfile",
            "/projects/ps-yeolab3/bay001/codebase/oligoCLIP/tests/ABC_multiplex/eclipse_multi.yaml",
            "-s",
            "/projects/ps-yeolab3/bay001/codebase/oligoCLIP/snakeCLIP.py",    
            "--use-singularity",

            "--directory",
            workdir,
        ], env=common.get_env())

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
