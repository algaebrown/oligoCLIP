import pandas as pd
from optparse import OptionParser
import warnings
import os
def parse_fastqc(filename):
    ''' A parser for fastQC .txt files
    input: str, full path for _fastqc.txt
    output:
        all_modules: a dictionary containing each QC metric's detail (in dataframe)
        psall: contain PASS/FAIL for each module (list in list)'''
    with open(filename) as f:
        module_lines = []
        psall = []
        all_modules = {}
        for line in f:
            
            if '>>END_MODULE' in line:
                all_modules[name]= module_lines
                psall.append([name, passfail])
                module_lines = []
                
            elif line.startswith('>>'): # start of new block
                
                name = line.split('\t')[0].replace('>>', '')
                passfail = line.split('\t')[1].rstrip()
                
            else:
                values = line.rstrip().split('\t')
                module_lines.append(values)
        
        all_modules[name]= module_lines
        psall.append([name, passfail])
        
        
        # make into dataframe
        for item in all_modules.keys():
            for i, line in enumerate(all_modules[item]):
                if '#' not in line[0]:
                    break # first line that is not column
            try:
                df = pd.DataFrame(all_modules[item][i:], columns = all_modules[item][i-1])
            except:
                df = pd.DataFrame(all_modules[item][i:])
            
            all_modules[item] = df
        
        return all_modules, psall

def option_parser():
    ''' return parser
    :return: OptionParser object
    '''
    usage = """
        this is a script to join many fastqc files into a csv
        python region_call.py -i FILE1.txt,FILE2.txt,FILE3.txt -b basic_stat.csv -p module.csv 
        """
    description = """FastqcIO Hsuan-lin Her 2020.
    
    example: python ~/projects/QC_tools/fastqc_io.py -i /home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_D2_IN.umi.r1.fqTrTr.sorted.fq.fastqc_data.txt,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_N2_IN.umi.r1.fqTrTr.sorted.fq.fastqc_data.txt,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_N1_IN.umi.r1.fqTrTr.sorted.fq.fastqc_data.txt,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_D1_IN.umi.r1.fqTrTr.sorted.fq.fastqc_data.txt,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_D2_IP.umi.r1.fqTrTr.sorted.fq.fastqc_data.txt,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_N2_IP.umi.r1.fqTrTr.sorted.fq.fastqc_data.txt,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_N1_IP.umi.r1.fqTrTr.sorted.fq.fastqc_data.txt,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_D1_IP.umi.r1.fqTrTr.sorted.fq.fastqc_data.txt -b basic.txt -p ps.txt
                   """
    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("--input", "-i", dest="input", help="comma-delimited filenames to _fastqc.txt", type="string")
    parser.add_option("--basic_statistics", "-b", dest="outfile", help="output file path to basic statistics, containing total number of reads, number that is poor quality and GC content", type="string")
    parser.add_option("--pass_fail_output", "-p",  dest="pass_fail_output", help="output file path to the .csv containing whether each QC metric failed or not", type="string")
    
    return parser

if __name__=='__main__':
    parser = option_parser()
    (options, args) = parser.parse_args()
    
    fastqc_files = options.input.split(' ')
    
    # check if file exist
    for f in fastqc_files:
        if not os.path.isfile(f):
            warnings.warn(f'{f} does not exist')
    
    # parse file
    all_data = [] # key -> value
    all_basic = []
    for file in fastqc_files:
        all_modules,pass_fail = parse_fastqc(file)
        df = pd.DataFrame(pass_fail).set_index(0)
        df.columns = [file]

        all_data.append(df)

        all_basic.append(all_modules['Basic Statistics'].set_index('#Measure')['Value'])
    
    # Join basic statistics
    all_data_basic_df=pd.concat(all_basic, axis = 1).T.reset_index()
    
    for col in ['Total Sequences','Sequences flagged as poor quality', '%GC']:
        all_data_basic_df[col] = all_data_basic_df[col].astype(float)
        print(f'======={col}==========')
        print(all_data_basic_df[col].describe())
    
    # for each module
    print('=== number of library failing each metric ===')
    all_data_df=pd.concat(all_data, axis = 1).T
    
    summary = all_data_df.apply(pd.Series.value_counts).fillna(0).T
    print(summary.loc[(summary['fail']>0)|(summary['warn']>0)])
    
    all_data_basic_df.to_csv(options.outfile)
    all_data_df.to_csv(options.pass_fail_output)
    
    
    



                