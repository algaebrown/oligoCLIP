import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_cutadapt_round_two():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/cutadapt_round_two/data")
        expected_path = PurePosixPath(".tests/unit/cutadapt_round_two/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("K562_rep1/fastqs/all.umi.fqTrTr.gz QC/K562_rep1.umi.r1.fqTrTr.metrics", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "K562_rep1/fastqs/all.umi.fqTrTr.gz",
            "QC/K562_rep1.umi.r1.fqTrTr.metrics",
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
        common.CutadaptOutputChecker(data_path, expected_path, workdir).check()
