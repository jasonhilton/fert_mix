#! /bin/bash

module load R/3.5.1

csv_file=$(ls run_specs/lko | tail -n 1)
runs=$(cat run_specs/lko/$csv_file | wc --lines)
arrayid=$(qsub PBS/run_stan_lko.pbs -v csv_file=$csv_file -t 1-$runs)
