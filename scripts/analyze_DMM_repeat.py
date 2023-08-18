from scipy.stats import dirichlet
import pandas as pd
from pathlib import Path
import numpy as np
import seaborn as sns
from scipy.stats import entropy
from scipy import stats
import matplotlib.pyplot as plt
import math
from scipy.stats import pearsonr
from scipy.spatial.distance import pdist, squareform
import sys
import tensorflow_probability as tfp
sns.set_palette('tab20c')
plt.style.use('seaborn-white')
    
def label_clusters_by_elbow(val, plot = False, ax = None, std = None):
    '''label clusters by model mean
    val: sorted model mean
    '''
    
    r_square = {}
    reg_params = {}
    for breakpoint in range(1, len(val)-1):
        
        left = val[:breakpoint+1]
        right = val[breakpoint:]
        
        left_x = np.arange(start = 0, stop = len(left), step = 1)
        left_reg = stats.linregress(left_x, left)
        right_x = np.arange(start = breakpoint, stop = breakpoint+len(right), step = 1)
        right_reg = stats.linregress(right_x, right)
        r_square[breakpoint] = left_reg.rvalue**2+right_reg.rvalue**2
        reg_params[breakpoint] = (left_x, left_reg, right_x, right_reg)
        
    elbow_point = max(r_square, key=r_square.get)
    left_x, left_reg, right_x, right_reg = reg_params[elbow_point]
    
    selected = val.iloc[elbow_point+1:]
    #print(val.loc[val-left_y>0].sort_values())
    
    
    if plot:
        if ax is None:
            f, ax = plt.subplots()
        left_x, left_reg, right_x, right_reg = reg_params[elbow_point]
        left_y = left_reg.intercept + left_reg.slope*left_x
        left_y_upper_bound = (left_reg.intercept+ left_reg.intercept_stderr) + (left_reg.slope+left_reg.stderr)*left_x
        left_y_lower_bound = (left_reg.intercept- left_reg.intercept_stderr) + (left_reg.slope-left_reg.stderr)*left_x
        ax.plot(left_x, left_y, 'tomato', label='left regression', lw = 3)
        ax.fill_between(left_x, left_y_upper_bound, left_y_lower_bound, color = 'pink', alpha = 0.5)
        ax.plot(right_x, right_reg.intercept + right_reg.slope*right_x, 'orchid', label='right regression', lw = 3)
        val.plot(marker = '+', color = 'grey', ax = ax, label = 'p_bar')
        if std is not None:
            ax.fill_between(np.arange(len(val)), val+3*std, val-3*std, color = 'skyblue', alpha = 0.5)
        ax.vlines(x = elbow_point, ymax = val.max(), ymin = 0, color = 'black', linestyle='dashed')
        ax.set_xticks(range(len(val)))
        ax.set_xticklabels(val.index, rotation = 90)
        ax.set_ylabel('$bar{p_{jk}}')
        ax.set_title(val.name)
    return selected

def annotate_clusters(model_mean, plot = False, model_std=None):
    cluster_assignment_df = []
    f, axes = plt.subplots(math.ceil(model_mean.shape[0]/2),2, figsize = (12,16))
    axes = axes.flatten()
    for rbp, ax in zip(model_mean.index, axes):

        val = model_mean.loc[rbp].sort_values()
        if model_std is not None:
            std = model_std.loc[rbp, val.index]
        else:
            std = None
        selected = label_clusters_by_elbow(val, plot = plot, ax = ax, std=std)
        selected = pd.Series([True]*len(selected), index = selected.index)
        selected.name = rbp
        cluster_assignment_df.append(selected)
    ax.legend()
    cluster_assignment_df = pd.concat(cluster_assignment_df, axis = 1).T
    cluster_assignment_df.fillna(False, inplace = True)
    plt.tight_layout()
    sns.despine()
    
    missing_annotation = list(set(model_mean.columns)-set(cluster_assignment_df.columns))
    cluster_assignment_df[missing_annotation]=False
    return cluster_assignment_df.T

def summaries_annotation_and_mean(anno, model_mean, model_var, read_dist):
    ''' create cluster summary statistics '''
    stats = []
    for index, row in anno.iterrows():
        #print(index)
        rbps = row.loc[row].index
        rbps_prefix_free = [r.split('.')[1] for r in rbps]
        mean = model_mean.loc[rbps, index].sum()
        var = model_var.loc[rbps, index].sum()
        bg_mean = read_dist.loc[rbps].sum()
        fc = mean/bg_mean
        
        comp_mean = model_mean[index].copy()
        comp_mean.index = [i.split('.')[1] for i in comp_mean.index]
        ri_contri = (comp_mean*np.log(comp_mean/read_dist))[rbps].sum()
        
        stat = [','.join(rbps_prefix_free), mean, var, bg_mean, fc, ri_contri]
        stats.append(stat)
    return pd.DataFrame(stats, index = anno.index, columns = ['RBP', 'RBP_mean', 'RBP_var', 'background_mean', 'fold_change', 'RI_contribution'])

def count_by_rbp(data_merged, anno, to_count = 'feature_type_top'):
    feature_counts = []
    for rbp in anno.columns:
        cnt = data_merged.loc[data_merged[rbp], to_count].value_counts()
        cnt.name = rbp.split('.')[1]
        feature_counts.append(cnt)
    feature_counts = pd.concat(feature_counts, axis = 1).fillna(0).T
    return feature_counts

def get_counts_by_rbp(data_merged, anno, data_col_to_merge = 'cluster', to_count = 'feature_type_top'):
    ''' count feature types by each RBP '''
    
    data_merged[anno.columns] = data_merged[anno.columns].fillna(False)
    
    feature_counts = count_by_rbp(data_merged, anno, to_count = to_count)
    return feature_counts

def compute_jaccard_index(identity_tbl):
    ''' use jaccard index to show how much binding site overlapped between RBPs'''
    d_condense = pdist(identity_tbl.T, 'jaccard')
    d = pd.DataFrame(1-squareform(d_condense), index = identity_tbl.columns, columns = identity_tbl.columns)
    return d


class Dirichlet_Mixture_Model:
    def __init__(self, alphas, weights, n):
        
        self.alphas = alphas # (alpha, beta) of each component, K*2 matrix
        self.weights = weights # vector of length K
        self.n = n
        
        self.distributions = []
        
        for a in self.alphas:
            self.distributions.append(tfp.distributions.DirichletMultinomial(n, a))
    def pmf(self,k):
        ''' return pmf of  mixture model'''
        p = 0
        for w, dist in zip(self.weights, self.distributions):
            p += w*dist.prob(k)
        return p
    
    def cdf(self,k):
        raise NotImplementedError
    def logpmf(self,k):
        return np.log(self.pmf(k))
    def pvalue(self,k):
        raise NotImplementedError

def DMM_bayes_factor(model_alphas, weights, rbps, components, raw_data, plot = True, nread = 30):
    '''dirichlet multinomial mixture log likelihood (BF) by aggregation property
    aggregate all the other proteins as noise
    calculate P(X|alt)/P(X|null)
    '''
    
    if weights.index[0] in model_alphas.columns:
        model_alphas.columns = [f'alpha.{c}' for c in model_alphas.columns]
    
    # marginalize into [rbp], [other]
    marginalized_alphas = pd.concat([model_alphas.loc[rbps].T, 
                                 model_alphas.loc[~model_alphas.index.isin(rbps)].sum(axis = 0)],
                               axis = 1)
    
    marginalized_counts = pd.concat([raw_data[rbps], 
                                 raw_data.loc[:, ~raw_data.columns.isin(rbps)].sum(axis = 1)],
                               axis = 1)
    
    
    comp_names = weights.loc[components].index # [V22]
    other_names = weights.loc[~weights.index.isin(components)].index # [all other comp]
    
    if 'alpha' in model_alphas.columns[0]:
        comp_names_alpha= [f'alpha.{c}' for c in comp_names]
        other_names_alpha = [f'alpha.{c}' for c in other_names]
    
    # reweight
    alt_w = weights.loc[comp_names, 'pi']
    alt_w = alt_w/alt_w.sum()
    alt = Dirichlet_Mixture_Model(marginalized_alphas.loc[comp_names_alpha].values,
                             alt_w,
                            n = marginalized_counts.sum(axis = 1).astype(float))
    
    
    null_w = weights.loc[other_names, 'pi']
    null_w = null_w/null_w.sum()
    null = Dirichlet_Mixture_Model(marginalized_alphas.loc[other_names_alpha].values,
                             null_w,
                             n = marginalized_counts.sum(axis = 1).astype(float))
    
    logL_comp = np.log(alt.pmf(marginalized_counts))
    logL_null = np.log(null.pmf(marginalized_counts))
    
    
    logLR = logL_comp-logL_null
    
    bf_df = pd.DataFrame([logL_null, logL_comp, logLR], 
                         index = ['logL_singlecomp', 'logL_comp', 'logLR'], 
                         columns = raw_data.index).T
    
    
    return bf_df
def filter_by_bf(bfs, individual_bfs, data, anno, comp_mapping = None, logbf_thres = 2):
    ''' filter hypothesis assignment by hypothesis-wise BF and individual BF 
    hypothesis BF: ex, How likely is this window being SF3B4+PRPF8 instead of all other hypothesis (including SF3B4 sole binding and PRPF8 sole binding)
    individual BF: ex, How likely is this window bound by SF3B4 (regardless of other partners' present) vs not being bound at all?
    '''
    data = data.copy()
    
    # filtering for combinatorial binding
    for index, row in data.iterrows():
        clus = row['cluster']
        if comp_mapping:
            try:
                evi = bfs.loc[index, comp_mapping[clus]]
                rbps = anno.loc[comp_mapping[clus]][anno.loc[comp_mapping[clus]]].index
            except:
                evi = 0
                rbps = None
        else:
            try:
                evi = bfs.loc[index, clus]
                rbps = anno.loc[clus][anno.loc[clus]].index
            except:
                evi = 0 # some of the clusters are not there anymore
                rbps = None
        
        if rbps is not None:
            if individual_bfs.loc[index, rbps].ge(logbf_thres).all():
                # if all RBPs are bound for specific component
                data.loc[index, 'logLR']=evi
            else:
                data.loc[index, 'logLR']=0
        else:
            data.loc[index, 'logLR']=0 # individual RBPs don't have evidence
        
        
        
    data.loc[data['logLR']>logbf_thres, 'BF_assignment'] = data.loc[data['logLR']>logbf_thres, 'cluster']
    
    return data



def mask_megaoutput(megaoutputs, mask):
    '''
    take megaoutput, modify logLR:{RBP} and binary {RBP} labels based on log_cdf and mask
    mask: True or false based on whether y_dev > zscore_cutoff * stdev
    '''
    # binary labels are gone
    megaoutputs_masked = megaoutputs.copy()
    common_index = list(set(megaoutputs_masked.index).intersection(set(mask.index)))
    megaoutputs_masked.loc[common_index, mask.columns] = megaoutputs.loc[common_index, mask.columns] & mask # being logLR > 2 and >2*stdev
    
    return megaoutputs_masked

if __name__ == '__main__':
    out_stem = sys.argv[1]
    basedir = Path(sys.argv[2])
    annotation_file = sys.argv[3]

    raw_counts_file = basedir / f'counts/repeats/megatables/name/{out_stem}.tsv.gz'
    annotation = pd.read_csv(annotation_file, sep = '\t', index_col = 0)
    name2family = annotation.set_index('repName')['repFamily'].to_dict()
    name2class = annotation.set_index('repName')['repClass'].to_dict()

    # constants
    logLR_threshold = 2
    ent_thres = 0.1
    fc_raw_thres = 1

    # fitted parameters and outputs
    data = pd.read_csv(basedir/'DMM_repeat/name'/f'{out_stem}.mixture_weight.tsv', sep = '\t', index_col = 0) # basedir/f'DMM/{out_stem}.mixture_weight.tsv'


    mixture_weight_only = data.loc[:, data.columns.str.startswith('X')]
    mixture_weight_only.columns
    data['cluster']=mixture_weight_only.idxmax(axis = 1)

    weights = pd.read_csv(basedir/'DMM_repeat/name'/f'{out_stem}.weights.tsv', sep = '\t', index_col = 0)
    weights.index = [f'X{i}' for i in weights.index]

    model_alphas = pd.read_csv(basedir/'DMM_repeat/name'/f'{out_stem}.alpha.tsv', sep = '\t',
                        index_col = 0) # RBP by components, B * K
    model_mean = model_alphas.div(model_alphas.sum(axis = 0), axis = 1)
    model_var = model_alphas.apply(lambda column: dirichlet(column).var(), axis = 0)
    model_std = np.sqrt(model_var)

    # raw_counts to calculate bayes factor
    raw_counts = pd.read_csv(raw_counts_file, sep = '\t', index_col = 0)
    raw_counts = raw_counts.loc[data.index] # filter for those that has been modelled (passing the total_read_threshold)
    mask = pd.read_csv(basedir / 'mask'/ f'{out_stem}.repeat_mask.csv', index_col = 0)

    # find total mapped reads as the null
    nread_per_window = raw_counts.sum(axis = 1)
    mapped_reads = raw_counts.sum(axis = 0)
    mapped_reads_fraction = mapped_reads.div(mapped_reads.sum())
    print('========Fraction mapped reads: \n ========',mapped_reads_fraction.sort_values())

    # visualization
    sns.clustermap(model_mean.T,
            cbar_kws = {'label': '\bar{p}'},
            metric = 'correlation', cmap = 'Greys', figsize = (4,4))
    plt.savefig(basedir /'DMM_repeat/name'/ f'{out_stem}.model_mean.pdf')


    # calculate FC over null for each component
    component_fc = (model_mean).div(mapped_reads_fraction, axis = 0).T
    # visualization
    sns.clustermap(component_fc,
            cbar_kws = {'label': 'FC over total mapped reads'},
            metric = 'correlation', cmap = 'Greys', figsize = (4,4))
    plt.savefig(basedir /'DMM_repeat/name'/ f'{out_stem}.component_fc.pdf')

    # annotate cluster: elbow method.
    anno = annotate_clusters(model_mean, plot = True, model_std = model_std)
    plt.savefig(basedir /'DMM_repeat/name'/ f'{out_stem}.elbow_labelling.pdf')

    # calculate entropy and filter, summarize clusters
    ent = model_mean.apply(lambda col: entropy(col, mapped_reads_fraction), axis = 0).sort_values()

    # contain how many masked regions
    data['contain_mask'] = data.index.isin(mask.index)
    fraction_mask = data.groupby(by = 'cluster')['contain_mask'].mean()

    # summarize clusters
    f, ax = plt.subplots(1,3, figsize = (9,3))
    annotation_summary = pd.concat([anno.sum(axis = 1), ent, fraction_mask], axis = 1)
    annotation_summary.columns = ['# RBP assigned', 'RI to total read distribution', 'fraction contain windows needing softmask']

    comp_stats = summaries_annotation_and_mean(anno, model_mean, model_var, mapped_reads_fraction)
    annotation_summary = comp_stats.merge(annotation_summary, left_index = True, right_index = True)
    annotation_summary['# RBP assigned'].sort_values().plot.barh(ax = ax[0], color = 'grey')

    annotation_summary['RI to total read distribution'].sort_values().plot.barh(ax = ax[1], color = 'grey')
    ax[1].set_ylabel('# relative entropy to\n  total reads distribution')
    sns.despine()

    anno.sum(axis = 0).sort_values().plot.barh(ax = ax[2], color = 'grey')
    ax[2].set_ylabel('# clusters assigned')
    plt.tight_layout()
    plt.savefig(basedir /'DMM_repeat/name'/ f'{out_stem}.cluster_summary.pdf')

    # low entropy is bad: should not filter because it still contains some individual binding sites
    # too_low_entropy = ent[ent<0.1].index.tolist()
    # annotation_summary['filtered']=annotation_summary.index.isin(too_low_entropy)
    # anno.loc[anno.index.isin(too_low_entropy)] = False

    annotation_summary.to_csv(basedir / f'{out_stem}.cluster_summary.csv')
    anno.to_csv(basedir /'DMM_repeat/name'/ f'{out_stem}.cluster_annotation_binary.csv')

    # calculate effect size: p_bar, fc_bar and p_bar's std
    ezik = data[model_mean.columns]
    p_bar = np.matmul(ezik, model_mean.T)
    p_bar.columns = model_mean.index
    p_raw = raw_counts.div(nread_per_window, axis = 0)
    fc_bar = p_bar.div(mapped_reads_fraction, axis = 1)
    fc_raw = p_raw.div(mapped_reads_fraction, axis = 1)

    var_bar = np.matmul(ezik, model_var.T)
    var_bar.columns = model_var.index
    std_bar = np.sqrt(var_bar)

    p_bar.to_csv(basedir /'DMM_repeat/name'/ f'{out_stem}.p_bar.csv')
    p_raw.to_csv(basedir /'DMM_repeat/name'/  f'{out_stem}.p_raw.csv')
    fc_bar.to_csv(basedir /'DMM_repeat/name'/  f'{out_stem}.fc_bar.csv')
    fc_raw.to_csv(basedir /'DMM_repeat/name'/  f'{out_stem}.fc_raw.csv')
    std_bar.to_csv(basedir /'DMM_repeat/name'/ f'{out_stem}.p_bar_std.csv')

    # plot model fit for p_bar and p_raw
    f, axes = plt.subplots(2, math.ceil(p_bar.shape[1]/2), figsize = (8, 5))
    axes = axes.flatten()

    for col, ax in zip(p_bar.columns, axes):
        ax.scatter(p_raw[col], p_bar.loc[p_raw.index, col], color = 'lightgrey', marker = '+')
        r, pval = pearsonr(p_raw[col], p_bar.loc[p_raw.index, col])
        ax.set_ylabel('p_bar_i')
        ax.set_xlabel('p_raw')
        ax.set_title(f'{col}\n r={r:.2f}\n p={pval:.2E}')
    sns.despine()
    plt.tight_layout()
    plt.savefig(basedir /'DMM_repeat/name'/ f'{out_stem}.p_fit.pdf')

    # annotate
    data['repFamily']=data.index.str.replace('_AS', '').map(name2family).tolist()
    data['repClass']=data.index.str.replace('_AS', '').map(name2class).tolist()
    # visualize regions distribution for each cluster
    col = 'repClass'
    cluster_count = data.pivot_table(index = 'cluster', columns = col, 
                                fill_value=0, aggfunc='size')
    cluster_frac = cluster_count.div(cluster_count.sum(axis = 1), axis = 0)

    # labelling 
    anno_as_color = anno.applymap(lambda ans: 'royalblue' if ans else 'white')
    anno_as_color.columns = [c.split('.')[1]+'(labels)' for c in anno_as_color.columns]
    sns.clustermap(cluster_frac*100,
            cmap = 'Greys', metric = 'cosine', cbar_kws ={'label': '%window', }, cbar_pos = (1,0.2,0.02,0.6),
            figsize = (5,3), xticklabels = 1, yticklabels = 1,
            row_colors = anno_as_color.T.sort_index().T)
    plt.savefig(basedir /'DMM_repeat/name'/ f'{out_stem}.cluster_repClass.pdf')
    cluster_count.to_csv(basedir /'DMM_repeat/name'/ f'{out_stem}.cluster_repClass.csv')

    col = 'repFamily'
    cluster_count = data.pivot_table(index = 'cluster', columns = col, 
                                fill_value=0, aggfunc='size')
    cluster_count.to_csv(basedir /'DMM_repeat/name'/ f'{out_stem}.cluster_{col}.csv')

    # Calculate Bayes Factor (hypothesis wise)
    comp_mapping = {}
    bfs_dmm = []
    for name, group in anno.groupby(by = list(anno.columns)):
        components = list(group.index)
        rbps = anno.columns[list(name)].tolist()
        if sum(name) > 0:
            print(components, rbps)
            bf_df = DMM_bayes_factor(model_alphas, weights, rbps, components, raw_counts)['logLR']
            bf_df.name = components[0]
            for comp in components:
                comp_mapping[comp]=components[0]
            bfs_dmm.append(bf_df)

    bfs_dmm = pd.concat(bfs_dmm, axis = 1)

    individual_bfs_dmm = []

    # test individual RBPs' binding
    for rbp in anno.columns:
        bound_components = anno[rbp][anno[rbp]].index
        individual_LR = DMM_bayes_factor(model_alphas, weights, [rbp], bound_components, raw_counts)['logLR']
        individual_LR.name = rbp
        individual_bfs_dmm.append(individual_LR)

    individual_bfs_dmm = pd.concat(individual_bfs_dmm, axis = 1)
    bfs_dmm.to_csv(basedir /'DMM_repeat/name'/ f'{out_stem}.BF_hypothesis.csv')
    individual_bfs_dmm.to_csv(basedir /'DMM_repeat/name'/ f'{out_stem}.BF_individual.csv') 

    # ====== By Hypothesis =====
    data = filter_by_bf(bfs_dmm, individual_bfs_dmm, 
                                        data = data, anno = anno, comp_mapping = comp_mapping)

    data_bf_dmm = data.merge(anno, left_on = 'BF_assignment', right_index = True, how = 'left')
    # rescue individual binding sites that were that detected by any hypothesis
    data_bf_dmm.loc[individual_bfs_dmm.index, anno.columns]=individual_bfs_dmm.ge(logLR_threshold)

    # mask it 
    data_bf_dmm_masked = mask_megaoutput(data_bf_dmm, mask)

    prefilter = get_counts_by_rbp(data_bf_dmm, anno, data_col_to_merge = 'BF_assignment',
                                                to_count = 'repClass'
                 )
    dmm_fcount_family = get_counts_by_rbp(data_bf_dmm_masked, anno, data_col_to_merge = 'BF_assignment',
                                                    to_count = 'repFamily')
    dmm_fcount_family.to_csv(basedir /'DMM_repeat/name'/ f'{out_stem}.rbp_family_type.csv')
    dmm_fcount_class = get_counts_by_rbp(data_bf_dmm_masked, anno, data_col_to_merge = 'BF_assignment',
                                                    to_count = 'repClass')
    dmm_fcount_family.to_csv(basedir /'DMM_repeat/name'/ f'{out_stem}.rbp_class_type.csv')

    diff = (prefilter-dmm_fcount_class).fillna(0)
    diff.loc[diff.sum(axis = 1).sort_values().index, diff.sum(axis = 0)>0].plot.bar(
    stacked = True, figsize = (3,3)
    )
    sns.despine()
    plt.savefig(basedir /'DMM_repeat/name'/ f'{out_stem}.n_masked.pdf')

    # counting
    dmm_fcount_family = get_counts_by_rbp(data_bf_dmm_masked, anno, data_col_to_merge = 'BF_assignment',
                                                    to_count = 'repFamily')
    dmm_fcount_family.to_csv(basedir /'DMM_repeat/name'/ f'{out_stem}.rbp_family_type.csv')
    dmm_fcount_class = get_counts_by_rbp(data_bf_dmm_masked, anno, data_col_to_merge = 'BF_assignment',
                                                    to_count = 'repClass')
    dmm_fcount_family.to_csv(basedir /'DMM_repeat/name'/ f'{out_stem}.rbp_class_type.csv')

    # calculate jaccard index
    if not data_bf_dmm_masked[anno.columns].sum(axis = 0).ge(1).all():
        # if not everything has at least 1 binding site
        cols_to_plot = anno.columns[data_bf_dmm_masked[anno.columns].sum(axis = 0).ge(1)]
    else:
        cols_to_plot = anno.columns


    if len(cols_to_plot)>1:
        d_dmm = compute_jaccard_index(data_bf_dmm_masked[anno.columns])
        cm=sns.clustermap(d_dmm, cmap = 'Greys', metric = 'correlation', figsize = (4,4),
                    vmax = 0.3)
        cm.cax.set_visible(False)
        plt.suptitle('ABC(Dirichlet Mixture Model: BF filtered)', y = 1)
        plt.savefig(basedir /'DMM_repeat/name'/ f'{out_stem}.jaccard_index.pdf')
    else:
        print('RBP with binding site:', cols_to_plot)

    try:
        dmm_fcount_family = dmm_fcount_family.loc[dmm_fcount_family.sum(axis = 1)>0, dmm_fcount_family.sum(axis = 0)>0]
        sns.clustermap(dmm_fcount_family.loc[:, ~(dmm_fcount_family.columns=='Simple_repeat')], 
            cmap = 'Greys', metric = 'correlation', figsize = (8,4), xticklabels = 1)
        plt.savefig(basedir /'DMM_repeat/name'/ f'{out_stem}.family_count.pdf')
        dmm_fcount_class = dmm_fcount_class.loc[dmm_fcount_class.sum(axis = 1)>0, dmm_fcount_class.sum(axis = 0)>0]
        sns.clustermap(dmm_fcount_class.loc[:, ~(dmm_fcount_class.columns=='Simple_repeat')], 
            cmap = 'Greys', metric = 'correlation', figsize = (5,4), xticklabels = 1)
        plt.savefig(basedir /'DMM_repeat/name'/ f'{out_stem}.class_count.pdf')
    except Exception as e:
        print(e)
    
    # ====== Generate Output ======
    # --- output by individual RBP ---
    columns = ['repFamily','repClass', 'logLR', 'cluster','BF_assignment']
    # output per RBP enrich windows
    data_bf_dmm_masked['name'] = data_bf_dmm_masked.index
    for c in anno.columns:

        enriched_windows = data_bf_dmm_masked.loc[data_bf_dmm_masked[c], columns]
        enriched_windows['p_bar']=p_bar[c]
        enriched_windows['p_raw']=p_raw[c]
        enriched_windows['fc_raw']=fc_raw[c]
        enriched_windows['fc_bar']=fc_bar[c]
        enriched_windows['p_bar_std']=std_bar[c]
        enriched_windows.to_csv(basedir / 'DMM_repeat/name'/f'{c}.enriched_windows.tsv', sep = '\t')
        print(f'found {c} enriched windows:', enriched_windows.shape)
        import os
        
    # --- output everything ---
    # output the full data
    p_bar.columns=[f'p_bar:{c}' for c in p_bar.columns]
    p_raw.columns=[f'p_raw:{c}' for c in p_raw.columns]
    fc_raw.columns=[f'fc_raw:{c}' for c in fc_raw.columns]
    fc_bar.columns=[f'fc_bar:{c}' for c in fc_bar.columns]
    std_bar.columns=[f'p_bar_std:{c}' for c in std_bar.columns]
    individual_bfs_dmm.columns=[f'logLR:{c}' for c in individual_bfs_dmm.columns]

    # concat everything and write together
    megaoutput = pd.concat([data_bf_dmm_masked, p_bar, p_raw, fc_raw, fc_bar, std_bar, individual_bfs_dmm], axis = 1)
    megaoutput.to_csv(basedir / 'DMM_repeat/name'/f'{out_stem}.megaoutputs.tsv', sep = '\t')

    data_bf_dmm.to_csv(basedir / 'DMM_repeat/name'/f'{out_stem}.megaoutputs_unmasked.tsv', sep = '\t')




















