import pandas as pd
from pybedtools import BedTool
import sys
def generate_variant_sequence(row):
    ref = row['REF']
    seq = row['seq']
    
    
    
    start = row['POS']-row['start']-1
    end = row['POS']-row['start']-1+len(ref)
    
    to_replace = seq[start:end]
    if to_replace != ref:
        print(to_replace, ref, start, end)
        return None
    new_seq = seq[:start]+row['ALT']+seq[end:]
    return new_seq

def reverse_complement(string):
    newstr = ''
    mapper={'A':'T',
            'T':'A',
            'C': 'G',
            'G':'C',
           '.': ''}
    for s in string[::-1]:
        newstr+=mapper[s]
    return newstr

def generate_variant_sequence_neg(row):
    ref = row['REF']
    seq = row['seq']
    
    
    
    start = row['end']-row['POS']
    end = row['end']-row['POS']-len(ref)
    
    to_replace = reverse_complement(seq[end+1:start+1])
    if to_replace != ref:
        print(to_replace, 'ref=',ref, start, end)
        return None
    new_seq = seq[:end+1]+reverse_complement(row['ALT'])+seq[start+1:]
    return new_seq

if __name__=='__main__':

    variants = pd.read_csv(sys.argv[1],
        sep = '\t',
        names = ['CHROM','POS','ID','REF', 'ALT', 'INFO/CLNDN', 'INFO/CLNVC', 'INFO/CLNSIG']
    )
    seq = pd.read_csv(sys.argv[2], sep = '\t',
                 names = ['chrom', 'name', 'seq', 'struct', 'label', 'start'])

    
    outf = sys.argv[4]

    if variants.empty:
        variant_df = pd.DataFrame(columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'name', 'variant_seq','feature_type_top', 'feature_types', 'gene_name',
        'transcript_types', 'transcript_type_top', 'INFO/CLNDN', 'INFO/CLNVC', 'INFO/CLNSIG'])
        variant_df.to_csv(outf)
    else:
        window_bed = BedTool(sys.argv[3])
        window_df = pd.read_csv(window_bed.fn, sep = '\t')

        # finding window name
        variants['POS-1']=variants['POS']-1

        df = BedTool.from_dataframe(variants[['CHROM','POS-1','POS','ID','REF', 'ALT','INFO/CLNDN', 'INFO/CLNVC', 'INFO/CLNSIG']]).intersect(
            window_bed, wb = True).to_dataframe(names = ['CHROM','POS-1', 'POS','ID','REF', 'ALT','INFO/CLNDN', 'INFO/CLNVC', 'INFO/CLNSIG']+window_df.columns.tolist())

        df['seq']=df['name'].map(seq.set_index('name')['seq'])
        df.dropna(subset = ['seq'],inplace = True)
        df = df.loc[~df['ALT'].str.contains('N')]
        df=df.loc[df['REF'].str.len()<10]

        pos = df.loc[df['strand']=='+']
        neg = df.loc[df['strand']=='-']
        
        # generate variant sequence
        # generate variant sequence
        if not pos.empty:
            pos['variant_seq']=pos.apply(generate_variant_sequence, axis = 1)
            
        else:
            pos['variant_seq']=None
        if not neg.empty:
            neg['variant_seq']=neg.apply(generate_variant_sequence_neg, axis = 1)
        else:
            neg['variant_seq']=None

        cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'name', 'variant_seq','feature_type_top', 'feature_types', 'gene_name',
        'transcript_types', 'transcript_type_top', 'INFO/CLNDN', 'INFO/CLNVC', 'INFO/CLNSIG']
        variant_df = pd.concat([pos[cols],
            neg[cols]],
            axis = 0)
        variant_df.dropna(subset = ['variant_seq'], inplace = True)

        variant_df.to_csv(outf)


