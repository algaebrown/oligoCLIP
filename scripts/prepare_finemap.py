import pyBigWig
import os
import sys
from optparse import OptionParser
from scipy.signal import find_peaks
from pybedtools import BedTool
import math
import pandas as pd
import numpy as np

# python scripts/region_breakdown.py \
# --ip1plus=/home/hsher/scratch/ENCODE_CITS/CITS/548_CLIP1.plus.bw \
# --ip1minus=/home/hsher/scratch/ENCODE_CITS/CITS/548_CLIP1.minus.bw \
# --ip2minus=/home/hsher/scratch/ENCODE_CITS/CITS/548_CLIP2.minus.bw \
# --ip2plus=/home/hsher/scratch/ENCODE_CITS/CITS/548_CLIP2.plus.bw \
# --inminus=/home/hsher/scratch/ENCODE_CITS/CITS/548_INPUT1.minus.bw \
# --inplus=/home/hsher/scratch/ENCODE_CITS/CITS/548_INPUT1.plus.bw \
# --region=/home/hsher/evan_regions/548_ZRANB2_K562.tsv \
# --bed=/home/hsher/548_ZRANB2_K562.bed
'''
python /home/hsher/projects/oligoCLIP/scripts/prepare_finemap.py \
        --ipminus /home/hsher/scratch/ABC_2rep_skipper_k562window/K562_rep4/bw/COV/RBFOX2.neg.bw \
        --ipplus /home/hsher/scratch/ABC_2rep_skipper_k562window/K562_rep4/bw/COV/RBFOX2.pos.bw \
        --inminus /home/hsher/scratch/ABC_2rep_skipper_k562window/K562_rep4/bw_bg/COV/RBFOX2.neg.bw \
        --inplus /home/hsher/scratch/ABC_2rep_skipper_k562window/K562_rep4/bw_bg/COV/RBFOX2.pos.bw \
        --region /home/hsher/scratch/ABC_2rep_skipper_k562window/internal_output/DMN/K562_rep4.RBFOX2.enriched_window.tsv \
        --bed test.bed

python ~/projects/oligoCLIP/scripts/prepare_finemap.py --ipminus K562_rep4/bw/COV/PRPF8.neg.bw --ipplus K562_rep4/bw/COV/PRPF8.pos.bw --inminus K562_rep4/bw_bg/COV/PRPF8.neg.bw --inplus K562_rep4/bw_bg/COV/PRPF8.pos.bw --region internal_output/DMN/K562_rep4.PRPF8.enriched_window.tsv --bed testprpf.bed
'''



def option_parser():
    ''' return parser
    :return: OptionParser object
    '''
    usage = """
        THIS IS CLIPPER FOR REGIONAL SHRINKING
        python shrink_region.py --ip <bamfile> --input <bamfile> --bed <bed> --out <fname>
        
        To shrink enriched regions into binding sites that we hope to contain motifs
Returns clusters (putative binding sites) and control clusters (size-matched, region-matched).
By comparing cluster to control clusters, we can calculated R values.
Clustering algorithm is DBSCAN
        """
    description = """CLIPper. Michael Lovci, Gabriel Pratt 2012, Hsuan-lin Her 2020.
                         CLIP peakfinder by shrinking enriched region."""
    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("--ipminus", dest="ipminus", help="CITS wigs", type="string")

    parser.add_option("--ipplus", dest="ipplus", help="CITS wigs", type="string")
    parser.add_option("--inminus", dest="inminus", help="CITS wigs", type="string")
    parser.add_option("--inplus", dest="inplus", help="CITS wigs", type="string")
    parser.add_option("--region", dest="regiontsv", help="Evan's region output", type="string")
    
    parser.add_option("--bed", "-b", dest="bed", help="bedfile output by region_call.py")

    return parser
    
    
class strand_specific_wig:
    def __init__(self, plus_wig, minus_wig):
        self.plus = pyBigWig.open(plus_wig)
        self.minus = pyBigWig.open(minus_wig)
    def fetch(self, chrom = None, start= None, end=None, strand= None, interval = None):
        ''' 
        return values for a bedtool interval or chrom, start, end, strand
        NOT 5' to 3'
        '''
        if interval:
            start = interval.start
            end = interval.end
            strand = interval.strand
            chrom = interval.chrom
        if strand == '-':
            icshape_data = self.minus
        else:
            icshape_data = self.plus
        values = np.array(icshape_data.values(chrom, start, end))
        if strand == '-':
            values = -values
        return values

def fetch_values(chrom, start, end, strand, IP_wig, IN_wig):
    
    # get crosslinking CITS
    IP_values = IP_wig.fetch(chrom = chrom, start = start, end = end, strand = strand)
    IN_values = IN_wig.fetch(chrom = chrom, start = start, end = end, strand = strand)
    
    # make into dataframe
    data = []
    starts = np.arange(start, end)
    ends = starts + 1

    data = pd.DataFrame([starts, ends, IP_values, IN_values], index = ['start', 'end', 'clip', 'input']).T
    data['chr'] = chrom
    data['strand'] = strand

    return data
    # col_names = c("chr","start","end","name","score","strand","window_n","input","clip")     



if __name__=='__main__':
    
    parser = option_parser()
    FDR_threshold = 0.2
    d_log_odds_threshold = 1
    
    (options, args) = parser.parse_args()

    print(options)

    # load evan's region
    df = pd.read_csv(options.regiontsv, sep = '\t')

    try:
        significant = df.loc[(df['qvalue'] < FDR_threshold) | (df['d_log_odds'] > d_log_odds_threshold)] 
    except:
        print('using all enriched windows')
        significant = df
    
    if significant.shape[0]==0:
        print('No significant enriched windows, no finemapping has to be done')
        with open(options.bed, 'w') as f:
            f.write('#zero enriched windows\n')
            f.write('#cols=chr,start,end,name,score,strand,window_n,input,clip\n')
    
    else:
        print('n_singificant_windows:', significant.shape[0])
        if 'chr' not in significant.columns:
            significant['chr'] = significant['chrom']
        

        if 'd_log_odds' in significant.columns:
            namecol = 'd_log_odds'
        elif 'mixture_model_score' in significant.columns:
            namecol = 'mixture_model_score'
        else:
            significant['qvalue1'] = significant['qvalue']
            namecol = 'qvalue1'
        # merge adjacent window
        cols = ['chr', 'start', 'end', 'name', 'qvalue', 'strand', namecol]
        sig_bed = BedTool.from_dataframe(significant[cols])
        merged_significant = sig_bed.sort().merge(s = True, c = [4,5,6,7], o = ['distinct', 'min', 'distinct', 'max'], d = 1
        ).to_dataframe(names = cols
        )

        print('n_merged_singificant_windows', merged_significant.shape)
        print(merged_significant.head())

        # get CITS
        IP_wig = strand_specific_wig(options.ipplus, 
                                options.ipminus)
        IN_wig = strand_specific_wig(options.inplus, 
                                options.inminus)
        
        # iterate over merged windows to fetch CITS/COV
        all_windows = []
        for index, row in merged_significant.iterrows():
            window_values = fetch_values(chrom = row['chr'], start = row['start'], end = row['end'], strand = row['strand'],
                IP_wig = IP_wig, IN_wig = IN_wig)
            # print(row)
            # print(window_values)
            window_values['window_n'] = row['name']
            window_values['score'] = row['qvalue']
            window_values['name'] = row[namecol]
            
            all_windows.append(window_values)

        all_window_values = pd.concat(all_windows, axis = 0).reset_index()

        all_window_values['start'] = all_window_values['start'].astype(int)
        all_window_values['end'] = all_window_values['end'].astype(int)
        all_window_values['input'].fillna(0, inplace = True)
        all_window_values['clip'].fillna(0, inplace = True)
        all_window_values['input'] = all_window_values['input'].astype(int)
        all_window_values['clip'] = all_window_values['clip'].astype(int)

        with open(options.bed, 'w') as f:
            f.write('#cols=chr,start,end,name,score,strand,window_n,input,clip\n')
            f.write(f'#name={namecol}\n')
            f.write(f'#score=qvalue\n')

            all_window_values[["chr","start","end","name","score","strand","window_n","input","clip"]].to_csv(f,
                sep = '\t', header = False, index = False)
        
        
        print(all_window_values.shape)
        print(all_window_values.head)