import sys
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import betabinom

class Beta_Mixture_Model:
    def __init__(self, alphas, betas, weights, n):
        self.alphas = alphas # (alpha, beta) of each component, K*2 matrix
        self.betas = betas
        self.weights = weights # vector of length K
        self.n = n
        
        self.distributions = []
        
        for a,b in zip(self.alphas, self.betas):
            self.distributions.append(betabinom(a = a, b = b, n = n))
        
    def pmf(self,k):
        ''' return pmf of  mixture model'''
        p = 0
        for w, dist in zip(self.weights, self.distributions):
            p += w*dist.pmf(k)
        return p
    
    def cdf(self,k):
        p = 0
        for w, dist in zip(self.weights, self.distributions):
            p += w*dist.cdf(k)
        return p
    def logpmf(self,k):
        return np.log(self.pmf(k))
    def pvalue(self,k):
        return 1-self.cdf(k)

if __name__=='__main__':
    basedir = Path(sys.argv[1]) # internal_output/DMN
    outstem = sys.argv[2] # K562_rep4.RBFOX2
    exp, rbp = outstem.split('.')
    raw_counts = pd.read_csv(sys.argv[3], sep = '\t') # basedir/f'internal_output/counts/genome/bgtables/internal/{pre}.{rbp}.tsv.gz
    annotation_df = pd.read_csv(sys.argv[4], sep = '\t')
    outdir = basedir

    # constants
    component_fold_threshold = 1
    fc_bar_threshold = 1
    logLR_threshold = 2
    FDR_cutoff = 0.2


    # find control columns
    col = raw_counts.columns
    control_col = [c for c in col if exp in c and c != outstem]
    assert len(control_col)==1
    control_col = control_col[0]

    # read files, weights of each component #pi
    component_weights = pd.read_csv(basedir/f'{outstem}.weights.tsv', sep = '\t',index_col = 0)
    component_weights.index = ['X'+str(i) for i in component_weights.index]

    # alpha 
    component_alpha = pd.read_csv(basedir/f'{outstem}.alpha.tsv', sep = '\t',index_col = 0)
    model_mean = component_alpha.div(component_alpha.sum(axis = 0), axis = 1)

    # E[z_ik]
    # weight of each window
    data = pd.read_csv(basedir/f'{outstem}.mixture_weight.tsv', sep = '\t',index_col = 0)
    data.rename({'Row.names':'name'}, inplace = True, axis = 1)
    to_rename = [c for c in data.columns if c.startswith('V')]
    new_name = [c.replace('V', 'X') for c in to_rename]
    data.rename(dict(zip(to_rename, new_name)), axis = 1, inplace = True)

    # find total mapped reads
    raw_counts = raw_counts.loc[raw_counts['name'].isin(data['name'])] # only those modelled
    raw_counts.set_index('name', inplace = True)
    counts = raw_counts[[outstem, control_col]]
    nread_per_window = counts.sum(axis = 1)
    mapped_reads = counts.sum(axis = 0)
    mapped_reads_fraction = mapped_reads.div(mapped_reads.sum())
    print('========Fraction mapped reads: \n ========',mapped_reads_fraction)

    if component_alpha.shape[1]>1:
        print('======== Multiple components detects, use mixture model to estimate binding========')
        # select components
        component_fc = model_mean.div(mapped_reads_fraction, axis = 0).T
        comp = component_fc.loc[component_fc[outstem]>component_fold_threshold].index
        not_bound_comp = component_fc.loc[component_fc[outstem]<=component_fold_threshold].index

        print(f'components with bindinig: {comp}')
        
        ezik = data.set_index('name')[model_mean.columns]
        ezik_marginalized = ezik[comp].sum(axis = 1)
        
        
        
        # p-value based on all non-bound components
        not_bound_mixture = Beta_Mixture_Model(component_alpha.loc[outstem, not_bound_comp],
                                            component_alpha.loc[control_col, not_bound_comp],
                                            component_weights.loc[not_bound_comp, 'pi'],
                                            n = nread_per_window)
        bound_mixture = Beta_Mixture_Model(component_alpha.loc[outstem, comp],
                                            component_alpha.loc[control_col, comp],
                                            component_weights.loc[comp, 'pi'],
                                            n = nread_per_window)
        p_alt = pd.Series(bound_mixture.logpmf(counts[outstem]), index = counts.index)
        p_null = pd.Series(not_bound_mixture.logpmf(counts[outstem]), index = counts.index)
        
        # binding score: estimated p from E[x_zik] and model_mean
        p_bar = np.matmul(ezik, model_mean.T)
        p_bar.columns = model_mean.index
        p_raw = counts[outstem]/nread_per_window
        
        
        
        p_df = pd.concat([p_raw, p_bar[outstem], ezik_marginalized, p_alt, p_null], axis = 1)
        p_df.columns = ['p_raw', 'p_bar', 'posterior_bound', 'log_L_bound', 'log_L_notbound']
        p_df['logLR']=p_df['log_L_bound']-p_df['log_L_notbound']
        
        f, ax=plt.subplots(4,1, figsize = (6,12))
        p_df.plot.scatter(x = 'p_raw', y = 'p_bar', color = 'grey', marker = '+', ax = ax[0])
        ax[0].set_title(outstem)
        
        
        p_df['fc_raw']=p_df['p_raw']/mapped_reads_fraction[outstem]
        p_df['fc_bar']=p_df['p_bar']/mapped_reads_fraction[outstem]
        
        p_df['fc_bar'].hist(bins = 100, color = 'grey', ax=ax[1])
        ax[1].set_xlabel('fc_bar')
        ax[1].set_ylabel('# windows')
        
        p_df.plot.scatter(x = 'posterior_bound', y = 'fc_bar', color = 'grey', marker = '+', ax = ax[2])
        
        
        
        p_df.plot.scatter(x = 'logLR', y = 'fc_bar', color = 'grey', marker = '+', ax = ax[3])
        plt.suptitle(outstem)
        plt.tight_layout()
        sns.despine()
        plt.savefig(outdir / f'{outstem}.model_output.pdf')
        
        plt.show()
        
        results = data.merge(p_df, left_on = 'name', right_index = True)
        results['enriched']=(results['fc_bar']>fc_bar_threshold)&(results['logLR']>logLR_threshold)
        
        # save metadata
        component_alpha.index = [f'alpha.{c}' for c in component_alpha.index]
        model_mean.index = [f'mean.{c}' for c in model_mean.index]
        component_fc.columns = [f'fc.{c}' for c in component_fc.columns]
        metadata = pd.concat([component_fc, model_mean.T, component_alpha.T], axis = 1)
        metadata['selected'] = metadata.index.isin(comp)
        metadata.to_csv(outdir / f'{outstem}.label_component.csv')
        
        bg_metadata = pd.concat([mapped_reads,mapped_reads_fraction], axis = 1)
        bg_metadata.columns = ['total_reads', 'fraction_reads']
        
        
        
    else:
        print(f'======== {outstem} has single component, fall back to hypothesis testing========')
        #data = calculate_pvalue(raw_counts, outstem, component_alpha.iloc[:, ], annotation_df)
        pvalue=1-betabinom.cdf(n = nread_per_window, a = component_alpha.loc[outstem].iloc[0], b = component_alpha.loc[control_col].iloc[0], k = counts[outstem])
        _, qvalue = fdrcorrection(pvalue)
        
        # fold change
        p_raw = raw_counts[outstem]/nread_per_window
        
        p_df = pd.DataFrame([pvalue, qvalue, p_raw], index = ['pvalue', 'qvalue', 'p_raw'], columns = counts.index).T
        p_df['fc_raw']=p_df['p_raw']/mapped_reads_fraction[outstem]
        
        results = data.merge(p_df, left_on = 'name', right_index = True)
        results['enriched']=(results['qvalue']<FDR_cutoff)

    enriched_windows = results.loc[results['enriched']]
    print(f'Finish testing, found enriched_windows: ', enriched_windows.shape[0])

    # save raw output
    results.to_csv(outdir / f'{outstem}.window_score.tsv', sep = '\t')
    enriched_windows.to_csv(outdir / f'{outstem}.enriched_windows.tsv', sep = '\t')

    # analysis
    fcount = results.groupby(by = 'enriched')['feature_type_top'].value_counts().unstack().fillna(0).T
    fcount['Positive rate'] = fcount[True]/(fcount[True]+fcount[False])
    fcount.to_csv(outdir / f'{outstem}.feature_type_summary.tsv', sep = '\t')

    gcount = results.groupby(by = 'enriched')['gene_type_top'].value_counts().unstack().fillna(0).T
    gcount['Positive rate'] = gcount[True]/(gcount[True]+gcount[False])
    gcount.to_csv(outdir / f'{outstem}.gene_type_summary.tsv', sep = '\t')

    tcount = results.groupby(by = 'enriched')['transcript_type_top'].value_counts().unstack().fillna(0).T
    tcount['Positive rate'] = tcount[True]/(tcount[True]+tcount[False])
    tcount.to_csv(outdir / f'{outstem}.transcript_type_summary.tsv', sep = '\t')

    f, ax=plt.subplots(3,2, figsize = (12,12))

    for i, cnt in enumerate([fcount, gcount, tcount]):
        cnt[True].sort_values().plot.barh(ax = ax[i,0], color = 'lightgrey')
        cnt['Positive rate'].sort_values().plot.barh(ax = ax[i,1], color = 'lightgrey')

    ax[i,0].set_xlabel('# enriched window')
    ax[i,1].set_xlabel('Positive rate')
    sns.despine()
    plt.tight_layout()
    plt.savefig(outdir / f'{outstem}.summaries.pdf')
    plt.show()

    # classification by feature
    binary_df = pd.DataFrame(False, index = results.index, columns = results['feature_type_top'].unique())
    for index, row in results.iterrows():
        binary_df.loc[index, row['feature_type_top'].split(':')] = True

    from sklearn.linear_model import LogisticRegression
    X = binary_df
    y = results.loc[binary_df.index, 'enriched']
    clf = LogisticRegression().fit(X, y)
    clf.score(X, y)

    pd.Series(clf.coef_[0], index = binary_df.columns).sort_values().plot.barh()
    plt.xlabel('Logistic Regression Coef')
    plt.tight_layout()
    plt.savefig(outdir / f'{outstem}.feature_logistic.pdf')
    plt.show()

    from sklearn.linear_model import RidgeClassifier
    X = binary_df
    y = results.loc[binary_df.index, 'enriched']
    clf = RidgeClassifier().fit(X, y)
    clf.score(X, y)

    pd.Series(clf.coef_[0], index = binary_df.columns).sort_values().plot.barh()
    plt.xlabel('Ridge Regression Coef')
    plt.tight_layout()
    plt.savefig(outdir / f'{outstem}.feature_ridge.pdf')
    plt.show()








