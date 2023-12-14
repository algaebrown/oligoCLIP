module load subreadfeaturecounts/1.5.3

gff=/tscc/projects/ps-yeolab4/software/skipper/1.0.0/bin/skipper/annotations/gencode.v38.annotation.k562.gt1.gff3.gz

full_gff=/tscc/projects/ps-yeolab4/software/skipper/1.0.0/bin/skipper/annotations/gencode.v38.annotation.gff3.gz

zcat $gff > ~/scratch/k562.gff
zcat $full_gff > ~/scratch/full.gff3

bam=~/scratch/oligo_PE_iter6_reseq/1021-Rep1/bams/CSTF2.rmDup.Aligned.sortedByCoord.out.bam
featureCounts -p -a ~/scratch/k562.gff -o ~/scratch/1021-Rep1.CSTF2.count $bam
featureCounts -p -a ~/scratch/full.gff3 -o ~/scratch/1021-Rep1.CSTF2.count.full $bam

bamr1=~/scratch/1021-Rep1.CSTF2.R1.bam
samtools view -hbf 64 $bam > $bamr1
featureCounts -p -a ~/scratch/k562.gff -o ~/scratch/1021-Rep1.CSTF2.r1.count $bamr1
featureCounts -p -a ~/scratch/full.gff3 -o ~/scratch/1021-Rep1.CSTF2.r1.count.full $bamr1


bam3=~/scratch/oligo_PE_iter5/oCLIP_25_1/bams/CSFT2.rmDup.Aligned.sortedByCoord.out.bam
featureCounts -p -a ~/scratch/k562.gff -o ~/scratch/iter5.CSFT2.count $bam3
featureCounts -p -a ~/scratch/full.gff3 -o ~/scratch/iter5.CSFT2.count.full $bam3

bamiter4=~/scratch/oligo_PE_iter4/Rep1/bams/PUM2.rmDup.Aligned.sortedByCoord.out.bam
featureCounts -p -a ~/scratch/k562.gff -o ~/scratch/iter4.PUM2.count $bamiter4
featureCounts -p -a ~/scratch/full.gff3 -o ~/scratch/iter4.PUM2.count.full $bamiter4

bamiter3=~/scratch/oligo_PE_iter3/oCLIP_Encode_Rep2_S2/bams/PUM2.rmDup.Aligned.sortedByCoord.out.bam
featureCounts -p -a ~/scratch/k562.gff -o ~/scratch/iter3.PUM2.count $bamiter3
featureCounts -p -a ~/scratch/full.gff3 -o ~/scratch/iter3.PUM2.count.full $bamiter3


bamspike=~/scratch/oligo_PE_iter7/1022-Rep1/bams/ctrlSpike.rmDup.Aligned.sortedByCoord.out.bam
featureCounts -p -a ~/scratch/k562.gff -o ~/scratch/spike.count $bamspike
featureCounts -p -a ~/scratch/full.gff3 -o ~/scratch/spike.count.full $bamspike

bam1plex=~/scratch/oligo_PE_iter7/1022-Rep1/bams/1plexRBFOX2.rmDup.Aligned.sortedByCoord.out.bam
featureCounts -p -a ~/scratch/k562.gff -o ~/scratch/1plexrbfox.count $bam1plex
featureCounts -p -a ~/scratch/full.gff3 -o ~/scratch/1plexrbfox.count.full $bam1plex

bam2plex=~/scratch/oligo_PE_iter7/1022-Rep1/bams/2plexRBFOX2.rmDup.Aligned.sortedByCoord.out.bam
featureCounts -p -a ~/scratch/k562.gff -o ~/scratch/2plexrbfox.count $bam2plex
featureCounts -p -a ~/scratch/full.gff3 -o ~/scratch/2plexrbfox.count.full $bam2plex

bambead=~/scratch/oligo_PE_iter7/1022-Rep1/bams/ctrlBead.rmDup.Aligned.sortedByCoord.out.bam
featureCounts -p -a ~/scratch/k562.gff -o ~/scratch/bead.count $bambead
featureCounts -p -a ~/scratch/full.gff3 -o ~/scratch/bead.count.full $bambead


bam2=~/scratch/ABC_2rep/K562_rep4/bams/PUM2.rmDup.Aligned.sortedByCoord.out.bam
featureCounts -a ~/scratch/k562.gff -o ~/scratch/K562_rep4.PUM2.count $bam2
featureCounts -a ~/scratch/full.gff3 -o ~/scratch/K562_rep4.PUM2.count.full $bam2




bam4=/tscc/nfs/home/hsher/scratch/katie_singleplex/output/bams/dedup/genome/RBFOX2_IP_1.genome.Aligned.sort.dedup.bam
featureCounts -a ~/scratch/k562.gff -o ~/scratch/katie_singleplex.RBFOX2_IP1 $bam4
featureCounts -a ~/scratch/full.gff3 -o ~/scratch/katie_singleplex.RBFOX2_IP1.full $bam4

bam5=/tscc/nfs/home/hsher/scratch/katie_singleplex/output/bams/dedup/genome/RBFOX2_IN_1.genome.Aligned.sort.dedup.bam
featureCounts -a ~/scratch/k562.gff -o ~/scratch/katie_singleplex.RBFOX2_IN1 $bam5
featureCounts -a ~/scratch/full.gff3 -o ~/scratch/katie_singleplex.RBFOX2_IN1.full $bam5

bam6=/oasis/tscc/scratch/eboyle/20220711_encode3_skipper/k562/output/bams/dedup/genome_R2/PUM1_IP_1.genome.Aligned.sort.dedup.R2.bam
featureCounts -a ~/scratch/k562.gff -o ~/scratch/eclip_pum1 $bam6
featureCounts -a ~/scratch/full.gff3 -o ~/scratch/eclip_pum1.full $bam6