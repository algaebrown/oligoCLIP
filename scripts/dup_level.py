import pandas as pd
import pysam
import sys

def get_mapped_read(bamfile):
    try:
        f = pysam.AlignmentFile(bamfile, "rb")
    except Exception as e:
        print(e)
        return 0
    return f.mapped
if __name__=='__main__':
    dup_files = sys.argv[1].split(' ')
    rmdup_files = sys.argv[2].split(' ')
    out = sys.argv[3]

    df = pd.DataFrame([dup_files, rmdup_files], index = ['dup_bam', 'rmdup_bam']).T

    df['before_dedup'] = df['dup_bam'].apply(get_mapped_read)
    df['after_dedup'] = df['rmdup_bam'].apply(get_mapped_read)

    df['percent unique fragment'] = df['after_dedup']/df['before_dedup']

    df.to_csv(out)