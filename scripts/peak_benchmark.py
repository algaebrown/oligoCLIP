# io
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
import os
import pandas as pd
import numpy as np

from pybedtools import BedTool

from plot_params import *
from color_mapper import *
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Set3.colors)
repeat_masker = BedTool('/home/hsher/gencode_coords/repeatmasker_ucsc.hg38.bed')

def read_annotated_peaks(path, filter_repeat = True, merge_bp_distance = 5, normalized = False):
    ''' read annotated peaks file as dataframe
    
    merges peaks within `merge_bp_distance`
    filter repeatmasker regions if `filter_repeat==True`
    return pd.DataFrame with columns -log10pval and score_bin
    '''
    normed_df = pd.read_csv(path, 
                 sep = '\t',
                names = ['chrom', 'start', 'end', 
                         'score', 'name', 'strand', 'gene_id',
                         'genename', 'region', 'detail'])
    
    normed_bed = BedTool.from_dataframe(normed_df).sort()
    
    if normalized:
        # the scores are logged
        normed_bed = normed_bed.merge(s = True, d = merge_bp_distance, c = [4,5,6,7,8,9,10]
                         , o = ['max', 'max', 'distinct', 'distinct', 'distinct', 'distinct', 'collapse'])
    else:
        # the scores are raw p-value
        if type(normed_df['score'][0])==str:
            
            normed_bed = normed_bed.merge(s = True, d = merge_bp_distance, c = [4,5,6,7,8,9,10]
                         , o = ['collapse', 'min', 'distinct', 'distinct', 'distinct', 'distinct', 'collapse'])
        else:
            normed_bed = normed_bed.merge(s = True, d = merge_bp_distance, c = [4,5,6,7,8,9,10]
                             , o = ['min', 'min', 'distinct', 'distinct', 'distinct', 'distinct', 'collapse'])

    
    if filter_repeat:
        before = normed_df.shape[0]
        normed_bed = normed_bed.intersect(repeat_masker, v = True, s = True)
        normed_df = normed_bed.to_dataframe(names = ['chrom', 'start', 'end', 'score', 'name', 'strand', 'gene_id',
                         'genename', 'region', 'detail'])
        after = normed_df.shape[0]
        
    
    
    if normalized:
        normed_df['-log10pval'] = normed_df['score']
    else:
        
        # can be either in name or score
        if normed_df['name'].max()<1:
            normed_df['-log10pval'] = -np.log10(normed_df['name']+normed_df.loc[normed_df['name']>0, 'name'].min())
        else:
            
            normed_df['-log10pval'] = -np.log10(normed_df['score']+normed_df.loc[normed_df['score']>0, 'score'].min())
    
    
    normed_df['score_bin'] = pd.cut(normed_df['-log10pval'], bins=np.arange(0,401),labels = False)
    
    return normed_df

###### metrics #########
def calculate_cdf(df, y = 'is_histone', x = 'score_bin'):
    ''' calcalate cdf of y variable, by gradually lowering x'''
    is_histone_cdf = np.cumsum(df.groupby(by = x)[y].sum()[::-1])
    total_peak_cdf = np.cumsum(df.groupby(by = x)[y].count()[::-1])
    
    return is_histone_cdf/total_peak_cdf

def get_prc_curve(y_test, y_score):
    lr_precision, lr_recall, thres =precision_recall_curve(y_test, y_score)
    lr_auc =auc(lr_recall, lr_precision)
    return lr_precision, lr_recall, lr_auc, thres

###### plotting ###########
def plot_CDF(normed_data_dict, unnormed_data_dict, dataset_to_plot, 
             color_mapper = color_mapper, 
             x = 'score_bin', y = 'is_histone',
            ylabel = 'fraction histone', xlabel = '-log10 pvalue'):
    ''' plotting CDF curve for normalize/unnormalized dataset based on y '''
    
    f, ax = plt.subplots(1,2,sharey = True)
    for label in dataset_to_plot:
        
        norm_cdf = calculate_cdf(normed_data_dict[label], y = y, x = x)
        unnorm_cdf = calculate_cdf(unnormed_data_dict[label], y = y, x = x)
        
        norm_cdf.plot(label = label, ax = ax[1], color = color_mapper(label))
        unnorm_cdf.plot(label = label, ax = ax[0], color = color_mapper(label))
        
        
    ax[1].legend(bbox_to_anchor = (1.5, 1))
    ax[0].set_ylabel(ylabel)
    ax[0].set_xlabel(xlabel)
    ax[1].set_xlabel(xlabel)
    ax[0].set_title('IP only')
    ax[1].set_title('normalized')
    sns.despine()
    
    
    
    
    
    
def plot_PRC(unnormed_data_dict, normed_data_dict, dataset_to_plot, 
             y_hat = '-log10pval', y_true = 'is_histone', color_mapper = color_mapper,
            region = None):
    ''' plot AUPRC curve '''

    f, ax = plt.subplots(1,2,sharey = True)
    for i in dataset_to_plot:

        unnormed = unnormed_data_dict[i]
        normed = normed_data_dict[i]
        color = color_mapper(i)
        
        if region is not None:
            unnormed = unnormed.loc[unnormed['region']==region]
            normed = normed.loc[normed['region']==region]

        try:
            pr, re, auc_, thres = get_prc_curve(unnormed[y_true], unnormed[y_hat])
            ax[0].plot(
            pr,
            re,
            color=color,
            label=f"(AUC = %0.2f)" % auc_,
            )
        except Exceptions as e:
            print(i, e)


        pr, re, auc_, thres = get_prc_curve(normed[y_true], normed[y_hat])


        ax[1].plot(
            pr,
            re,
            color=color,
            label=f"{i} (AUC = {auc_:2f})",
        )


    ax[0].set_ylabel('Recall')
    ax[0].set_xlabel('Precision')
    ax[1].set_xlabel('Precision')
    ax[1].legend(bbox_to_anchor = (1.2, 0.8))
    ax[0].set_title('IP only')
    ax[1].set_title('normalize to IN')
    plt.suptitle(f'AUPRC (positive = {y_true}; y_hat={y_hat})')
    sns.despine()