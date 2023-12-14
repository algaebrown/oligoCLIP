cd ~/scratch/oligo_PE_iter6_reseq/

# original comment
bedtools bamtobed -i 1021-Rep2/bams/PRPF8.rmDup.Aligned.sortedByCoord.out.bam \
    | awk '($1 != "chrEBV") && ($4 !~ "/2$")' \
    | bedtools flank -s -l 1 -r 0 -g /tscc/nfs/home/hsher/gencode_coords/GRCh38.primary_assembly.chrom.sizes -i - \
    | bedtools shift -p -1 -m 1 -g /tscc/nfs/home/hsher/gencode_coords/GRCh38.primary_assembly.chrom.sizes -i - \
    | bedtools coverage -counts -s -a /tscc/projects/ps-yeolab4/software/skipper/1.0.0/bin/skipper/annotations/gencode.v38.annotation.k562_totalrna.gt1.tiled_partition.bed.gz -b - \
    | cut -f 7 \
    | awk 'BEGIN {print "1021-Rep2.PRPF8"} {print}' > counts/genome/vectors/1021-Rep2.PRPF8.counts;

bedtools bamtobed -i 1021-Rep2/bams/PRPF8.rmDup.Aligned.sortedByCoord.out.bam \
    | awk '($1 != "chrEBV") && ($4 !~ "/2$")' \
    | bedtools flank -s -l 1 -r 0 -g /tscc/nfs/home/hsher/gencode_coords/GRCh38.primary_assembly.chrom.sizes -i - \
    | bedtools shift -p -1 -m 1 -g /tscc/nfs/home/hsher/gencode_coords/GRCh38.primary_assembly.chrom.sizes -i - \
    > ~/scratch/shifted.bed 

genome_window=/tscc/projects/ps-yeolab4/software/skipper/1.0.0/bin/skipper/annotations/gencode.v38.annotation.k562_totalrna.gt1.tiled_partition.bed.gz 
repeat_window=/tscc/projects/ps-yeolab3/eboyle/encode/pipeline/05_20220720/repeatmasker.grch38.sort.unique.bed.gz

bedtools intersect -s -v -a ~/scratch/shifted.bed \
-b $genome_window > ~/scratch/notgenome.bed

bedtools intersect -s -v -a ~/scratch/notgenome.bed \
-b $repeat_window > ~/scratch/aintnothing.bed

#wc -l 
751131 /tscc/nfs/home/hsher/scratch/shifted.bed
662429 /tscc/nfs/home/hsher/scratch/notgenome.bed
340085 /tscc/nfs/home/hsher/scratch/aintnothing.bed

# MAPQ
cut -d $'\t' -f 5 /tscc/nfs/home/hsher/scratch/aintnothing.bed | sort | uniq -c
313826 0 # multimapping reads
   8160 1
  12002 255
   6097 3

bedtools makewindows -g /tscc/nfs/home/hsher/gencode_coords/GRCh38.primary_assembly.chrom.sizes -w 500 > ~/scratch/genome500.bed
# 
bedtools coverage -counts -a ~/scratch/genome500.bed -b /tscc/nfs/home/hsher/scratch/aintnothing.bed > ~/scratch/genome500.count

awk ' $4 > 1000 ' ~/scratch/genome500.count > ~/scratch/highcovshit.bed

cut -d $'\t' -f 4  ~/scratch/highcovshit.bed | paste -sd+ | bc
#290499 explains (85%)

cut -d $'\t' -f 1  ~/scratch/highcovshit.bed | sort | uniq -c
    #   3 chr21
    #  31 GL000220.1
    #  32 KI270733.1