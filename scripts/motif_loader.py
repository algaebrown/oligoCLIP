# loads RBNS and SELEX motifs given RBP name

import pickle
from scipy.stats import entropy
import pandas as pd
import numpy as np
import os
rbns_status = pd.read_csv('/tscc/nfs/home/hsher/projects/rbns_thermo/RBNS_status.csv', names = ['rbp', 'status'])
rbns_status['status'] = rbns_status['status'].fillna('NAN').str.upper()
rbns_file_dict = pickle.load(open('/tscc/nfs/home/hsher/projects/rbns_thermo/RBNS_filename.pickle', 'rb'))
rbp_with_rbns = rbns_status.loc[rbns_status['status']=='PASS', 'rbp'].tolist()

selex_dir = '/tscc/nfs/home/hsher/selex_motif'


def get_RBNS_motif_PWM(rbp):
    ''' return PWM matrix '''
    motif_fname = rbns_file_dict[rbp]
    if len(motif_fname)>1:
        fname = [f for f in motif_fname if '0nM' not in f][0]
    else:
        fname = motif_fname[0]

    motif = pd.read_csv(fname, sep = '\t', skiprows = 1, index_col = 0)
    
    return motif

def get_SELEX_motif_PWM(rbp):
    
    linear_motifs = os.listdir(os.path.join(selex_dir, 'linear'))
    motif_path = [os.path.join(selex_dir, 'linear', f) for f in linear_motifs if rbp in f]
    if len(motif_path) == 0:
        struct_motifs = os.listdir(os.path.join(selex_dir, 'struct'))
        print(struct_motifs)
        motif_path = [os.path.join(selex_dir, 'struct', f) for f in struct_motifs if rbp in f]
        if len(motif_path)==0:
            return None
    
    motif = pd.read_csv(motif_path[0], index_col = 0)
    
    return motif

class Motif:
    def __init__(self, pwm, pseudocount = 0.000001):
        self.pwm = pwm
        
        self.filter_by_entropy()
        self.pwm += pseudocount

        # clean up index, make U to T (same as bedtools)
        print(self.pwm)
        self.pwm.columns = ['T' if i == 'U' else i for i in list(self.pwm.columns)]
        
    def filter_by_entropy(self, thres = 1):
        '''calculate entropy for each position '''
        m_entropy = np.apply_along_axis(arr = self.pwm.values, func1d = entropy, axis = 1)
        self.pwm = self.pwm.loc[m_entropy <thres]
        self.pwm.reset_index(inplace = True)
    def score_samelen_seq(self, seq):
        ''' return log likelihood score of a sequence '''
        score = 0
        assert len(seq) == self.pwm.shape[0]
        for index, s in enumerate(seq):
            score += np.log(self.pwm.loc[index, s]) # the higher the better (higher prob)
        return score
    def score_entire_seq(self, seq):
        ''' return the scores of every kmer in sequence '''
        k = self.pwm.shape[0]
        scores = []
        for i in range(len(seq)-k+1):
            substr = seq[i:i+k]
            substr_score = self.score_samelen_seq(substr)
            scores.append(substr_score)
        return scores
    def max_scoring_substring(self, seq):
        scores = self.score_entire_seq(seq)
        max_score = max(scores)
        index_start = scores.index(max_score)
        index_end = index_start+ self.pwm.shape[0]
        return max_score, index_start, index_end