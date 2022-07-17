# summarize peak.annotate.py
# how many significant peaks
# how many RNA species
# what types of regions

from peak_benchmark import read_annotated_peaks
from optparse import OptionParser
import logging
import pandas as pd
import os

def optparser():
    parser = OptionParser()
    parser.add_option("-b", "--bed", dest="annotated_bed",
                    help="annotated bed file from annotator")
    parser.add_option('-n', '--is_normalized', default = False, action = "store_true")
    parser.add_option("-f", "--filter_repeat",
                    action="store_true", dest="filter_repeat", default=False,
                    help="filter_out_repeat_masker")
    parser.add_option("-m", "--merge_bp_distance",
                    default=0,
                    help="merge peak within how many bp")
    parser.add_option("-t", "--log10pval_thres",
                    default=3,
                    help="-log10 pvalue threshold to filter peaks")
    parser.add_option("--l2fold_change_threshold",
                    default=3,
                    help="log2 FC threshold to filter peaks. Does not apply for unnormalized peaks")
    parser.add_option('-o', '--outfile', dest = "outfile", help = 'output path of summary statistics')
    (options, args) = parser.parse_args()
    return options, args

def is_histone(string):
    ''' based on genename, determine whether it is a histone gene'''
    if 'HIST' in string or string.startswith('H1') | string.startswith('H2') | string.startswith('H3') | string.startswith('H4'):
        return True
    else:
        return False


if __name__ == '__main__':
    options, args = optparser()
    logging.info(options)
    peak_df = read_annotated_peaks(options.annotated_bed,
                                   filter_repeat = options.filter_repeat,
                                   merge_bp_distance = options.merge_bp_distance,
                                   normalized = options.is_normalized
                                    )
    n_prefilter = peak_df.shape[0]
    logging.info('pre-filter:',n_prefilter)
    # filter pvalue
    if options.is_normalized:
        # has both pvalue and fc
        peak_df = peak_df.loc[(peak_df['name']>options.l2fold_change_threshold)&
                                (peak_df['-log10pval']>options.log10pval_thres)
                                ]
    else:
        peak_df = peak_df.loc[(peak_df['-log10pval']>options.log10pval_thres)
                                ]
    n_postfilter = peak_df.shape[0]
    logging.info('post-filter:',n_postfilter)

    region_count = peak_df['region'].value_counts()
    region_perc = region_count.div(region_count.sum())
    region_perc.index = [r+'_fraction' for r in region_perc.index]

    transcript_type = peak_df['detail'].str.split(':', expand = True)[8].value_counts()
    transcript_type_perc = transcript_type.div(transcript_type.sum())
    transcript_type_perc.index = [r+'_fraction' for r in transcript_type_perc.index]

    n_snoRNA = peak_df['genename'].str.startswith('SNOR').sum()
    n_mtRNA = peak_df['genename'].str.startswith('MT-').sum()
    n_histone = peak_df['genename'].apply(is_histone).sum()

    summary = pd.concat([region_count, region_perc, transcript_type, transcript_type_perc], axis = 0)
    summary['before_filter']=n_prefilter
    summary['after_filter'] = n_postfilter
    summary['n_snoRNA'] = n_snoRNA
    summary['n_mtRNA'] = n_mtRNA
    summary['n_histone'] = n_histone

    with open(options.outfile, 'w') as f:
        f.write('#'+str(options)+'\n')
        f.write(os.path.basename(options.annotated_bed))
        summary.to_csv(f)






    