rule get_gnomAD_snp:


    shell:
        """
        bcftools query -R /home/hsher/scratch/ABC_2rep/DMM/finemapping/mapped_sites/CITS/K562_rep6.RBFOX2.finemapped_windows.bed.gz -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz
        """

rule fetch_sequence:

    