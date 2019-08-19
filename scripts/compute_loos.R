library(rstan)
library(purrr)
library(dplyr)
library(loo)
library(magrittr)
library(tidyr)
library(dplyr)
library(stringr)
source("R/assessment_functions.R")
source("R/convergence_functions.R")
source("R/loo_lko_functions.R")


cmd_arg <- commandArgs(trailingOnly = T)
run_dir <- cmd_arg[1]

if (is.na(run_dir)){
  # the latest one
  run_dir <- tail(sort(list.files("mcmc_results")),1)
  run_dir <- paste0("mcmc_results/", run_dir)
}

time_stamp <- stringi::stri_extract(run_dir, regex="20[0-9]{6}_[0-9]{6}")
countries <- list.files(run_dir, pattern="[A-Z]{3,}")

# cycle through all models, computing and saving out loos into the processed
# directory
walk(countries, compute_country_loos, run_dir)

## Create a list of runs to refit the model with points with high pareto ks 
# left out.
run_spec <- map_df(countries, get_country_lko_run_list, run_dir)

run_spec %<>% mutate(base_config="lko_hb") %>% select(-file_name)
dir.create(file.path("run_specs", "lko"), recursive = T)
readr::write_csv(run_spec, 
                 path=file.path("run_specs", "lko",
                                paste0("run_",
                                       time_stamp,
                                       ".csv")))