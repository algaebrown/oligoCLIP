
window=/tscc/projects/ps-yeolab4/software/skipper/1.0.0/bin/skipper/annotations/gencode.v38.annotation.k562_totalrna.gt1.tiled_partition.bed.gz
window_peak_format=/tscc/nfs/home/hsher/scratch/k562_window.peaks

# make peak file: peakID, chr, start, end, strand
zcat $window |  awk -F $'\t' '{print $4, $1,$2,$3,$6}' OFS='\t' > $window_peak_format

module load homer
pum2_motif=/tscc/nfs/home/hsher/scratch/ABC_2rep/beta-mixture_CC/homer/finemapped_results/COV/K562_rep6.PUM2/homerMotifs.motifs8
pum2_motif_dir=/tscc/nfs/home/hsher/scratch/ABC_2rep/beta-mixture_CC/homer/finemapped_results/COV/K562_rep6.PUM2/
rbfox2_motif=/tscc/nfs/home/hsher/scratch/ABC_2rep/beta-mixture_CC/homer/finemapped_results/COV/HEK293_rep1.RBFOX2/homerMotifs.motifs6
rbfox2_motif_dir=/tscc/nfs/home/hsher/scratch/ABC_2rep/beta-mixture_CC/homer/finemapped_results/COV/


findMotifsGenome.pl $window_peak_format hg38 $pum2_motif_dir -find $pum2_motif -size given > /tscc/nfs/home/hsher/scratch/k562_window_pum2.homer.motif
findMotifsGenome.pl $window_peak_format hg38 $rbfox2_motif_dir -find $rbfox2_motif -size given > /tscc/nfs/home/hsher/scratch/k562_window_rbfox2.homer.motif