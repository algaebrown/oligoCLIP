rep=oligo_PE_iter6

rm ~/scratch/$rep/*/bams -rf
rm ~/scratch/$rep/QC/mapping*
rm ~/scratch/$rep/QC/summary.csv
rm ~/scratch/$rep/QC/read_count -rf


# change chromosome assembly in config file

snakemake -s snakeOligoCLIP_PE.smk -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" \
--configfile config/preprocess_config/oligope_iter6.yaml --use-conda --conda-prefix /tscc/nfs/home/hsher/snakeconda QC/summary.csv -n

snakemake -s snakeABC_SE.smk -j 12 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo -e {params.error_out_file} -o {params.out_file}" \
--configfile config/preprocess_config/oligose_k562_noalt.yaml --use-conda --conda-prefix /tscc/nfs/home/hsher/snakeconda QC/summary.csv -n