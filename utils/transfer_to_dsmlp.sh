rsync -r ~/scratch/ABC_2rep/variants/CITS/K562_rep6.*.csv hsher@dsmlp-login.ucsd.edu:bin/PrismNet/data/variants -v
rsync -r ~/scratch/ABC_2rep/variants_clinvar/CITS/K562_rep6.*.csv hsher@dsmlp-login.ucsd.edu:bin/PrismNet/data/variant_clinvar -v

# sync back inference results
rsync -rv hsher@dsmlp-login.ucsd.edu:bin/PrismNet/data/variant_clinvar/*variant_score.csv ~/scratch/variant_clinvar_score
rsync -rv hsher@dsmlp-login.ucsd.edu:bin/PrismNet/data/variants/*variant_score.csv ~/scratch/variant_score
