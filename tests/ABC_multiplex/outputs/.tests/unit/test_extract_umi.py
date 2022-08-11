import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_extract_umi():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/extract_umi/data")
        expected_path = PurePosixPath(".tests/unit/extract_umi/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("K562_rep1/fastqs/all.umi.fq.gz QC/K562_rep1.all.umi.metrics", file=sys.stderr)
        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "K562_rep1/fastqs/all.umi.fq.gz",
            "QC/K562_rep1.all.umi.metrics",
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
        common.FastqOutputChecker(data_path, expected_path, workdir).check()
