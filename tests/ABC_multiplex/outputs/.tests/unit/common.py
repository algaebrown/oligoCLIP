"""
Common code for unit testing of rules generated with Snakemake 7.3.8.
"""

from pathlib import Path
import subprocess as sp
import os

from parse_cutadapt import parse_new_cutadapt_file

def get_env():
    # use docker login to get around any bandwidth limitations that might arise.
    env = dict(os.environ, SINGULARITY_BINDPATH='/home,/projects,/oasis') # , SINGULARITY_DOCKER_USERNAME=os.environ['SINGULARITY_DOCKER_USERNAME'], SINGULARITY_DOCKER_PASSWORD=os.environ['SINGULARITY_DOCKER_PASSWORD'])
    
    return env
    
class OutputChecker:
    def __init__(self, data_path, expected_path, workdir):
        self.data_path = data_path
        self.expected_path = expected_path
        self.workdir = workdir

    def ignore_extensions(self, f):
        '''
        Unit tests check outputs against expected. Currently snakeCLIP does not 
        list these intermediates, so unit tests don't expect these files. 
        Benchmarked files also get erroneously thrown into here, which tests 
        do not create.
        '''
        return False
    
    def check(self):
        input_files = set(
            (Path(path) / f).relative_to(self.data_path)
            for path, subdirs, files in os.walk(self.data_path)
            for f in files
        )
        expected_files = set(
            (Path(path) / f).relative_to(self.expected_path)
            for path, subdirs, files in os.walk(self.expected_path)
            for f in files
        )
        unexpected_files = set()
        for path, subdirs, files in os.walk(self.workdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.workdir)
                if str(f).startswith(".snakemake"):
                    continue
                if f in expected_files:
                    self.compare_files(self.workdir / f, self.expected_path / f)
                elif f in input_files:
                    # ignore input files
                    pass
                elif self.ignore_extensions(f):
                    pass
                else:
                    unexpected_files.add(f)
        if unexpected_files:
            raise ValueError(
                "Unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files)))
                )
            )
        print(expected_files)
        print(unexpected_files)
    def compare_files(self, generated_file, expected_file):
        sp.check_output(["cmp", generated_file, expected_file])


class IndexOutputChecker(OutputChecker):

    def ignore_extensions(self, f):
        for ext in ['.index_bam.txt']:
            if str(f).endswith(ext):
                return True
        return False

    
class FastqcOutputChecker(OutputChecker):
    
    def ignore_extensions(self, f):
        for ext in ['.fo', '.zip', '.png', 'summary.txt', 'fastqc_report.html']:
            if str(f).endswith(ext):
                return True
        if '.cache/' in str(f) or 'benchmarks/' in str(f) or '.java/' in str(f):
            return True
        return False
    
    def compare_files(self, generated_file, expected_file):
        if str(generated_file).endswith('.txt'):  # ignores checking html file
            sp.check_output(["cmp", generated_file, expected_file])
        return 0

    
class DedupOutputChecker(OutputChecker):

    def ignore_extensions(self, f):
        for ext in ['_edit_distance.tsv', '_per_umi.tsv', '_per_umi_per_position.tsv']:
            if str(f).endswith(ext):
                return True
        if '.cache/' in str(f) or '/benchmarks/' in str(f):
            return True
        return False
        
    
class FastqOutputChecker(OutputChecker):
    def ignore_extensions(self, f):
        if '.cache/' in str(f) or 'benchmarks/' in str(f):
            return True
        return False
    
    def compare_files(self, generated_file, expected_file):
        if str(generated_file).endswith('.gz'):
            generated_fqlen = int(sp.check_output([f'zcat {generated_file} | wc -l'], shell=True).rstrip())
            expected_fqlen = int(sp.check_output([f'zcat {generated_file} | wc -l'], shell=True).rstrip())
            assert generated_fqlen == expected_fqlen
        elif str(generated_file).endswith('.fastq') or str(generated_file).endswith('.fq'):
            generated_fqlen = int(sp.check_output(['wc', '-l', generated_file]).rstrip())
            expected_fqlen = int(sp.check_output(['wc', '-l', expected_file]).rstrip())
            assert generated_fqlen == expected_fqlen
            
        return 0
    
    
class CutadaptOutputChecker(OutputChecker):

    def ignore_extensions(self, f):
        for ext in ['.txt']:
            if str(f).endswith(ext):
                return True
        return False
    
    def compare_files(self, generated_file, expected_file):
        if str(generated_file).endswith('.metrics'):
            generated_metrics = parse_new_cutadapt_file(generated_file)
            expected_metrics = parse_new_cutadapt_file(expected_file)
            print(generated_metrics)
            print(expected_metrics)
            assert generated_metrics['Reads Written'] == expected_metrics['Reads Written']
        elif str(generated_file).endswith('.gz'):
            generated_fqlen = int(sp.check_output([f'zcat {generated_file} | wc -l'], shell=True).rstrip())
            expected_fqlen = int(sp.check_output([f'zcat {expected_file} | wc -l'], shell=True).rstrip())
            assert generated_fqlen == expected_fqlen
        return 0


class SnakeleafOutputChecker(OutputChecker):
    def compare_files(self, generated_file, expected_file):
        # Snakeleaf just has date and time, no use comparing to expected
        # which will have a different datetime anyway...
        return 0


class BamOutputChecker(OutputChecker):
    
    def ignore_extensions(self, f):
        '''
        Unit tests check outputs against expected. Currently snakeCLIP does not 
        list these intermediates, so unit tests don't expect these files. 
        Benchmarked files also get erroneously thrown into here, which tests 
        do not create.
        '''
        for ext in ['.Log.out', '.Log.progress.out', '.SJ.out.tab', 
            '.align_reads.txt', '.align_dros_reads.txt'
        ]:
            if str(f).endswith(ext):
                return True
        return False

    def compare_files(self, generated_file, expected_file):
        mapped_reads = {}
        
        if str(generated_file).endswith('.bam'):
            for f in [generated_file, expected_file]:
                mapped_reads[f] = int(sp.check_output([
                    '/projects/ps-yeolab4/software/yeolabconda3/envs/samtools-1.15/bin/samtools', 
                    'view', '-c', '-F', '4', f'{f}'
                ]).rstrip())
            
            assert mapped_reads[generated_file] == mapped_reads[expected_file]
        return 0