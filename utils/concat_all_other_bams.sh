module load samtools
OUTPUT=/tscc/nfs/home/hsher/scratch/ABC_2rep_K562_rep6.bam
SORTED_OUTPUT=/tscc/nfs/home/hsher/scratch/ABC_2rep_K562_rep6.sorted.bam
samtools merge $OUTPUT /tscc/nfs/home/hsher/scratch/ABC_2rep/K562_rep6/bams/*rmDup*.out.bam
samtools sort $OUTPUT | samtools view -Sb - > $SORTED_OUTPUT
samtools index $SORTED_OUTPUT