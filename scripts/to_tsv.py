import pandas as pd
import sys
if __name__=="__main__":
    in_csv = sys.argv[1]
    out_tsv = sys.argv[2]
    df = pd.read_csv(in_csv)

    df[['rbp', 'barcode']].to_csv(out_tsv, sep = '\t',
                                index=False,
                                header = False)