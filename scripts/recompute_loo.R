library(rstan)
library(purrr)
library(dplyr)
library(loo)
library(stringr)
library(magrittr)
library(tidyr)
library(dplyr)
library(matrixStats)

source("R/loo_lko_functions.R")

cmd_arg <- commandArgs(trailingOnly = T)
run_dir <- cmd_arg[1]


if (is.na(run_dir)){
  # the latest one
  run_dir <- tail(sort(list.files("mcmc_results/lko")),1)
  run_dir <- paste0("mcmc_results/lko/", run_dir)
}
run_time <- stringi::stri_extract(run_dir, regex="20[0-9]{6}_[0-9]{6}")


run_spec_lko <- read.csv(paste0("run_specs/lko/run_", run_time, ".csv")) %>%
  as_tibble() %>% mutate(Index = 1:n())



countries <- list.files(run_dir, pattern="[A-Z]{3,}")

loos <- list()

assess <-readRDS("processed/assess_out.Rds")

#runs <- run_spec_lko %>% filter(country==cntry)
for (i in 1:4){
  cntry <- countries[i]
  message("Country = ", cntry)
  # assess <- readRDS(file.path("processed", paste0("assess_", cntry, ".rds"))) %>%
  #   filter(dispersion=="variable")
  
  assess_cntry <- assess %>% filter(country==cntry)
  loo_cntry <- map_df(seq_len(dim(assess_cntry)[1]),
      function(j){
        get_correct_loo(assess_cntry[j,],run_spec_lko, run_dir)
      }
  )
  loos[[i]] <- loo_cntry
}

loo_df <- bind_rows(loos)  
loo_df %<>% mutate(Components=paste0(str_to_title(param_type_1),
                                     " / ",
                                     str_to_title(param_type_2)))

saveRDS(loo_df, "processed/corrected_loos.rds")




