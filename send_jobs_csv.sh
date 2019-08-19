#! /bin/bash

module load R/3.5.1

Rscript scripts/make_csv_batch.R base_hb

csv_file=$(ls run_specs | tail -n 1)
runs=$(cat run_specs/$csv_file | wc --lines)
arrayid=$(qsub PBS/run_stan_csv.pbs -v csv_file=$csv_file -t 1-$runs)

qsub -W depend=afterokarray:$arrayid PBS/run_loos.pbs