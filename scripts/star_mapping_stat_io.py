import pandas as pd
from optparse import OptionParser
import warnings
import os

mapping_df = []
def read_mapstats(fname):
    '''  '''
    with open(fname) as f:
        lines = f.readlines()
    df = pd.DataFrame([l.rstrip().split('|\t') for l in lines ]).set_index(0)
    df.index = [i.replace('  ', '') for i in df.index]
    return df
def concat_mapstats(filenames, names = None):
    ''' given a list of filenames, read all of them and concat into a pd.DataFrame '''
    mapping_df = []
    for i,f in enumerate(filenames):
        try:
            df = read_mapstats(f)
            df.loc['STAR Log filename'] = [f]
            df.columns = [os.path.basename(f).split('.')[1]]
            if names:
                df.loc['name'] = names[i]
            mapping_df.append(df)
        except Exception as e:
            print(e)
            print(f)
    total = pd.concat(mapping_df, axis = 1).T
    
    total.columns = [c.rstrip(' ').strip(' ') for c in total.columns]
    
    return change_data_type(total)

def change_data_type(total):
    ''' given a dataframe, change each columns to appropriate datatype '''
    # cleaning dataframe
    for col in total.columns:
        if 'Number' or 'Average' in col:
            try:
                total[col] = total[col].astype(float)
            except Exception as e:
                print(col, e)
        if '%' in col:
            total[col] = total[col].str.replace('%', '').astype(float)
    return total
    
def option_parser():
    ''' return parser
    :return: OptionParser object
    '''
    usage = """
        this is a script to join many fastqc files into a csv
        python star_mapping_stat_io.py -i FILE1.txt,FILE2.txt,FILE3.txt -o stat.csv
        """
    description = """FastqcIO Hsuan-lin Her 2020.
    
    example:  python ~/projects/QC_tools/star_mapping_stat_io.py -i /home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_N2_IP.umi.r1.fq.repeat-unmapped.sorted.STARLog.final.out,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_D2_IN.umi.r1.fq.repeat-unmapped.sorted.STARLog.final.out,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_D2_IP.umi.r1.fq.repeat-unmapped.sorted.STARLog.final.out,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_N2_IN.umi.r1.fq.repeat-unmapped.sorted.STARLog.final.out,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_N1_IP.umi.r1.fq.repeat-unmapped.sorted.STARLog.final.out,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_D1_IN.umi.r1.fq.repeat-unmapped.sorted.STARLog.final.out,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_D1_IP.umi.r1.fq.repeat-unmapped.sorted.STARLog.final.out,/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_N1_IN.umi.r1.fq.repeat-unmapped.sorted.STARLog.final.out -o STAR_stas.csv
                   """
    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("--input", "-i", dest="input", help="comma-delimited filenames to .Log.final.out from STAR aligner", type="string")
    parser.add_option("--basic_statistics", "-o", dest="outfile", help="output to the concat csv", type="string")
    
    parser.add_option("--min_read", dest="nread", help="threshold of # of read that will produce warning", type="int", default=10**7)
    parser.add_option("--unique_map_rate",dest="unique", help="threshold of unique mapping rate that will produce warning", type="int", default=35)
    parser.add_option("--min_read_length", dest="rlen", help="threshold of read length that will trigger warning", type="int", default=55)
    
    
    return parser



class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    
if __name__=='__main__':
    parser = option_parser()
    (options, args) = parser.parse_args()
    
    log_files = options.input.split(' ')
    
    min_reads=options.nread
    unique_map_rate=options.unique
    min_read_length=options.rlen
    
    total = concat_mapstats(log_files)
    
    # check input reads
    if total['Number of input reads'].le(min_reads).any():
        print(f"{bcolors.WARNING}===Some have less than {min_reads} reads ==={bcolors.ENDC}")
        
        print(total.loc[total['Number of input reads'].le(min_reads), 'Number of input reads'])
    
    if total['Uniquely mapped reads %'].le(unique_map_rate).any():
        print(f"{bcolors.WARNING}===Some have low unique mapping rate ==={bcolors.ENDC}")
        
        print(total.loc[total['Uniquely mapped reads %'].le(unique_map_rate), 'Uniquely mapped reads %'])
        
        print(f'{bcolors.HEADER}Here are the reasons why they are not mapping well{bcolors.ENDC}')
        unmap_reasons = ['% of reads mapped to too many loci', 
                    '% of reads unmapped: too many mismatches', 
                    '% of reads unmapped: too short',
                   '% of reads unmapped: other']
        print(total.loc[total['Uniquely mapped reads %'].le(unique_map_rate), unmap_reasons])
        
        print('''if you see TOO SHORT for most of them, it means alignment is too short. 
              You might not be trimming adaptors correctly. check fastQC adaptor content.
              Also check the unmapped reads and see if they come from bacteria.....
              Also check the input read length. ''')
        
    if total['Average input read length'].le(min_read_length).any():
        print(f"{bcolors.WARNING}===Some read lengths are very short. Are you trimming them badly? ==={bcolors.ENDC}")
        print("check fastQC outputs to see how long they originally are")
        print("check cutadapt command")
        print("check RNase condition, cleanup steps and tapestation")
        
        print(total.loc[total['Average input read length'].le(min_read_length)].index.tolist())
        
    total.to_csv(options.outfile)
        
    
    
    
    

