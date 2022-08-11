import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_align_to_GENOME():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/align_to_GENOME/data")
        expected_path = PurePosixPath(".tests/unit/align_to_GENOME/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("K562_rep1/bams/genome/IGF2BP2.genome-mapped.Aligned.out.bam K562_rep1/bams/genome/IGF2BP2.genome-mapped.Unmapped.out.mate1 K562_rep1/bams/genome/IGF2BP2.genome-mapped.Log.final.out", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "K562_rep1/bams/genome/IGF2BP2.genome-mapped.Aligned.out.bam",
            "K562_rep1/bams/genome/IGF2BP2.genome-mapped.Unmapped.out.mate1",
            "K562_rep1/bams/genome/IGF2BP2.genome-mapped.Log.final.out",
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
        common.BamOutputChecker(data_path, expected_path, workdir).check()
