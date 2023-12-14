'''
python /tscc/nfs/home/hsher/projects/oligoCLIP/scripts/analyze_betabinom_mixture.py \
    /tscc/nfs/home/hsher/scratch/ABC_2rep_skipper_k562window/internal_output/DMN \
    K562_rep4.RBFOX2 \
    /tscc/nfs/home/hsher/scratch/ABC_2rep_skipper_k562window/internal_output/counts/genome/bgtables/internal/K562.RBFOX2.tsv.gz \
    /tscc/projects/ps-yeolab4/software/skipper/1.0.0/bin/skipper/annotations/gencode.v38.annotation.k562_totalrna.gt1.tiled_partition.features.tsv.gz

python /tscc/nfs/home/hsher/projects/oligoCLIP/scripts/analyze_betabinom_mixture.py \
    /tscc/nfs/home/hsher/scratch/ABC_2rep_skipper_k562window/internal_output/DMN \
    K562_rep4.PUM2 \
    /tscc/nfs/home/hsher/scratch/ABC_2rep_skipper_k562window/internal_output/counts/genome/bgtables/internal/K562.PUM2.tsv.gz \
    /tscc/projects/ps-yeolab4/software/skipper/1.0.0/bin/skipper/annotations/gencode.v38.annotation.k562_totalrna.gt1.tiled_partition.features.tsv.gz


'''
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import betabinom


def label_component(model_alphas, null_alphas, label, min_fold = 3):
    ''' label which components are RBP-enriched by the alpha(generalization of p)
    model_alphas: alphas from beta-binom mixture model.
    null_alphas: alphas from single component beta-binom model
    label: string to search for enrichment,
    min_fold: f, for a component to be included, its RI has to be > f*(min(all RI))
    '''
    
    # p*log(p/q)
    model_alphas_ri = model_alphas.apply(
        lambda col: (col/np.sum(col))*np.log2((col/np.sum(col))/(null_alphas.iloc[:,0]/null_alphas.iloc[:,0].sum())), axis = 0)
        
    ri_cutoff = model_alphas_ri.sum(axis = 0).min()*min_fold
    return model_alphas_ri, model_alphas_ri.loc[:,(model_alphas_ri.idxmax(axis = 0).str.contains(label))&(model_alphas_ri.sum(axis = 0)>ri_cutoff)].columns.tolist()

def beta_binom_cdf(counts, alphas, column):
    n = counts.sum()
    dist = betabinom(n, alphas[0], alphas[1])
    
    return dist.cdf(counts[column])

def calculate_pvalue(raw_counts, col_to_search, null_alpha, annotation_df):
    ''' calculate p-value using beta-binom mixture given alphas '''
    data_col = null_alpha.index.tolist()

    # find unique count combinations
    rcount = raw_counts[data_col].drop_duplicates()
    rcount['pvalue']=rcount.apply(
        lambda row:1-beta_binom_cdf(row, 
                                  null_alpha.loc[data_col].tolist(), 
                                  column = col_to_search)
                , axis = 1)
    
    # map p-value back to 
    raw_counts = raw_counts.merge(rcount, left_on = data_col,
                                  right_on = data_col,
                             how = 'left')

    mcol = ['chr', 'start', 'end', 'name', 'strand']
    mcol1 = ['chrom', 'start', 'end', 'name', 'strand']
    raw_counts = raw_counts.merge(annotation_df, left_on = mcol, right_on = mcol1)
    _, raw_counts['qvalue'] = fdrcorrection(raw_counts['pvalue'])
    
    return raw_counts
if __name__=='__main__':
    basedir = Path(sys.argv[1]) # internal_output/DMN
    outstem = sys.argv[2] # K562_rep4.RBFOX2
    raw_counts = pd.read_csv(sys.argv[3], sep = '\t') # basedir/f'internal_output/counts/genome/bgtables/internal/{pre}.{rbp}.tsv.gz
    annotation_df = pd.read_csv(sys.argv[4], sep = '\t')

    
    
    score_cutoff = 0.6
    fdr_cutoff = 0.2
    nread_thres = 10

    # read files, weights of each component
    component_weights = pd.read_csv(basedir/f'{outstem}.weights.tsv', sep = '\t',index_col = 0)
    component_weights.index = ['X'+str(i) for i in component_weights.index]
    
    # alphas
    null_alpha = pd.read_csv(basedir/f'{outstem}.null.alpha.tsv', sep = '\t',index_col = 0)
    component_alpha = pd.read_csv(basedir/f'{outstem}.alpha.tsv', sep = '\t',index_col = 0)
    
    # weight of each window
    data = pd.read_csv(basedir/f'{outstem}.mixture_weight.tsv', sep = '\t',index_col = 0)
    data.rename({'Row.names':'name'}, inplace = True, axis = 1)
    to_rename = [c for c in data.columns if c.startswith('V')]
    new_name = [c.replace('V', 'X') for c in to_rename]
    data.rename(dict(zip(to_rename, new_name)), axis = 1, inplace = True)
    
    # label components
    model_alphas_ri, comp =  label_component(component_alpha, null_alpha, outstem)

    # filter raw_counts
    raw_counts = raw_counts.loc[raw_counts[null_alpha.index].sum(axis = 1)> nread_thres]

    try:
        sns.clustermap(model_alphas_ri.T,cbar_kws = {'label': 'p*log(p/q)'},
                      metric = 'correlation', cmap = 'Greys', figsize = (4,4))
        plt.savefig(basedir / f'{outstem}.model_ri.pdf')
    
    except Exception as e:
        print(e)

    ########## WHEN THERE IS MORE THAN 1 COMP, USE MIXTURE WEIGHT ##############
    if component_alpha.shape[1]>1:
        # output which component is selected and why
        metadata = component_alpha.copy()
        metadata.index = [c+'_alpha' for c in metadata.index]

        model_alphas_ri.index = [c+'_ri' for c in metadata.index]
        model_alphas_ri.loc['RI'] = model_alphas_ri.sum(axis = 0)

        metadata = pd.concat([metadata.T, model_alphas_ri.T, component_weights], axis = 1)
        metadata.loc[comp, 'selected'] = True
        metadata['selected'].fillna(False, inplace = True)
        metadata.to_csv(basedir / f'{outstem}.label_component.csv')

        # calculate qvalue from the most null component
        raw_counts = calculate_pvalue(raw_counts, outstem, component_alpha[model_alphas_ri.sum(axis = 0).idxmin()], annotation_df)
        
        # label
        data['mixture_model_score'] = data[comp].mean(axis = 1)
        data['pvalue'] = data['name'].map(raw_counts.set_index('name')['pvalue'])
        data['qvalue'] = data['name'].map(raw_counts.set_index('name')['qvalue'])
        if len(comp)>0:
            enriched_features = data.loc[data[comp].ge(score_cutoff).any(axis = 1)]
        else:
            enriched_features = data.loc[data['qvalue']<fdr_cutoff]

        
    else:
        print(f'{outstem} has single component, fall back to hypothesis testing')
        data = calculate_pvalue(raw_counts, outstem, null_alpha['single_component_weight'], annotation_df)
        enriched_features = data.loc[data['qvalue']<fdr_cutoff]
        

    # save raw output
    data.to_csv(basedir / f'{outstem}.window_score.tsv', sep = '\t')
    enriched_features.to_csv(basedir / f'{outstem}.enriched_window.tsv', sep = '\t')


