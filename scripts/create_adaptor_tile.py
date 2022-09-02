from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import sys

def create_tiling_sequence(ada, outf, length=10):
    records = []
    for i in range(len(ada)-length):
        
        record = SeqRecord(
        Seq(ada[i:i+length]),
        id=f"{i}",
        )
        
        records.append(record)
    SeqIO.write(records, outf, format = 'fasta')

if __name__=='__main__':
    ada = sys.argv[1]
    outf = sys.argv[2]
    k = int(sys.argv[3])

    create_tiling_sequence(ada, outf, length = k)