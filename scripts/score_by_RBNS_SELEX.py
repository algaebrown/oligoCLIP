from motif_loader import *
import sys
import pandas as pd
from Bio import SeqIO

if __name__ == '__main__':
    rbns_selex_name = sys.argv[1]
    fa = sys.argv[2]
    out = sys.argv[3]

    try:
        pwm = get_RBNS_motif_PWM(rbns_selex_name)
    except:
        pwm = get_SELEX_motif_PWM(rbns_selex_name)
    
    motif = Motif(pwm)

    # read all the sequence and score them
    data = []
    for record in SeqIO.parse(fa, "fasta"):
        try:
            lls, rela_start, rela_end = motif.max_scoring_substring(str(record.seq).upper())
            
            data.append([lls, rela_start, rela_end, record.id, record.seq])
        except Exception as e:
            print(len(record.seq))
    
    df = pd.DataFrame(data, columns = ['motif_LLS', 'motif_start', 'motif_end', 'peak_id', 'seq'])
    df.to_csv(out)