# Configuration file 
import pandas as pd


#======= exp 1 ===========
fastqs=['/tscc/projects/ps-yeolab5/hsher/ABC009_1_S1_R1_001.fastq.gz',
'/tscc/projects/ps-yeolab5/hsher/ABC009_2_S2_R1_001.fastq.gz',
'/tscc/projects/ps-yeolab5/hsher/ABC009_3_S3_R1_001.fastq.gz',
'/tscc/projects/ps-yeolab5/hsher/ABC009_4_S4_R1_001.fastq.gz',
'/tscc/projects/ps-yeolab5/hsher/ABC011_5_S5_R1_001.trimmed.fastq.gz',
'/tscc/projects/ps-yeolab5/hsher/ABC011_6_S6_R1_001.trimmed.fastq.gz',
'/tscc/projects/ps-yeolab5/hsher/ABC011_7_S7_R1_001.trimmed.fastq.gz',
'/tscc/projects/ps-yeolab5/hsher/ABC011_8_S8_R1_001.trimmed.fastq.gz',
'/tscc/projects/ps-yeolab5/hsher/ABC011_1_S1_R1_001.trimmed.fastq.gz',
'/tscc/projects/ps-yeolab5/hsher/ABC011_2_S2_R1_001.trimmed.fastq.gz',
'/tscc/projects/ps-yeolab5/hsher/ABC011_3_S3_R1_001.trimmed.fastq.gz',
'/tscc/projects/ps-yeolab5/hsher/ABC011_4_S4_R1_001.trimmed.fastq.gz']

names = ['K562_rep1',
'K562_rep2',
'HEK293_rep1',
'HEK293_rep2',
'K562_SLBP_rep1',
'K562_SLBP_rep2',
'K562_RBFOX2_rep1',
'K562_RBFOX2_rep2',
'K562_rep3',
'K562_rep4',
'K562_rep5',
'K562_rep6']



fastq_menifest = pd.DataFrame([fastqs, names], index = ['fastq', 'libname']).T

print(fastq_menifest)


fastq_menifest.iloc[:4].to_csv('multiplex1.csv')
fastq_menifest.iloc[4:6].to_csv('singleplex_slbp.csv')
fastq_menifest.iloc[6:8].to_csv('singleplex_rbfox.csv')
fastq_menifest.iloc[8:].to_csv('multiplex2.csv')