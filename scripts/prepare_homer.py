import pandas as pd
import sys


if __name__=='__main__':
    enriched_window = sys.argv[1]
    tested_window = sys.argv[2]
    d_log_odds_cutoff = float(sys.argv[3])
    annotation = sys.argv[4]
    output_prefix = sys.argv[5]

    enriched_df = pd.read_csv(enriched_window, sep = '\t')
    tested_df = pd.read_csv(tested_window, sep = '\t')
    anno_df = pd.read_csv(annotation, sep = '\t')

    # map features to annotation df
    tested_df['feature_type_top'] = tested_df['name'].map(anno_df.set_index('row_id')['feature_type_top'])

    # filter enriched
    enriched_df = enriched_df.loc[enriched_df['d_log_odds']> d_log_odds_cutoff]
    
    # remove enriched from tested
    tested_df = tested_df.loc[~tested_df['name'].isin(enriched_df['name'].tolist())]

    assert len(set(enriched_df['name']).intersection(set(tested_df['name'])))==0

    # homer columns
    homer_cols = ['name', 'chr', 'start', 'end', 'strand']
    
    # group by region
    for region in ['UTR3', 'UTR5', 'CDS', 'INTRON']:
        feat_enrich_sub = enriched_df.loc[enriched_df['feature_type_top'].str.contains(region)]
        feat_test_sub = tested_df.loc[tested_df['feature_type_top'].str.contains(region)]

        # make format homer likes
        feat_enrich_sub[homer_cols].to_csv(f'{output_prefix}.{region}.enrich.homer', sep = '\t', header = None, index = None)
        feat_test_sub[homer_cols].to_csv(f'{output_prefix}.{region}.tested.homer', sep = '\t', header = None, index = None)
    
    enriched_df[homer_cols].to_csv(f'{output_prefix}.ALL.enrich.homer', sep = '\t', header = None, index = None)
    tested_df[homer_cols].to_csv(f'{output_prefix}.ALL.tested.homer', sep = '\t', header = None, index = None)