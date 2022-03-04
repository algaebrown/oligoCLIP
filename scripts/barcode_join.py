import pandas as pd
import sys

if __name__=='__main__':
    barcode_csv = sys.argv[1]
    whitelist = sys.argv[2]
    out = sys.argv[3]

    barcode=pd.read_csv(barcode_csv, header = None, sep = '\t', names = ['RBP', 'barcode'])
    whitelist=pd.read_csv(whitelist, sep = '\t', comment = '#', 
        header = None, 
        names = ['barcode', 'alternative', 'occurence', 'alt_occurence'],
        index_col = None)
    
    merged = barcode.merge(whitelist, left_on = 'barcode', right_on = 'barcode')

    merged.to_csv(out)
