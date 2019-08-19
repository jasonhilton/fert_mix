library(rstan)
library(dplyr)
library(loo)
library(magrittr)
library(tidyr)
library(dplyr)
library(purrr)
source("R/rate_reconstruction.R")
source("R/assessment_functions.R")
source("R/convergence_functions.R")
source("R/data_functions.R")

cmd_arg <- commandArgs(trailingOnly = T)
run_dir <- cmd_arg[1]

if (is.na(run_dir)){
  # the latest one
  run_dir <- tail(sort(list.files("mcmc_results")),1)
  run_dir <- paste0("mcmc_results/", run_dir)
}


# Collate all assessment files -------------------------------------------------
ass_files <- list.files( run_dir, pattern="assessment")

assessment_df  <- map_df(ass_files, function(fl){
  
  readRDS(file.path(run_dir, fl))
})

assessment_df<-assessment_df[!duplicated(assessment_df),]

# Find multimodel models -------------------------------------------------------
nonconverged_rows <- assessment_df %>% filter(!converged_rhat) %>% 
  filter(dispersion=="variable")


get_model_obj <- function(nonconverged_row, run_dir){
  non_con_model <- readRDS(file.path(run_dir,
                                     nonconverged_row$country,
                                     paste("mixture",
                                           nonconverged_row$param_type_1,
                                           nonconverged_row$param_type_2,
                                           nonconverged_row$dispersion,
                                           "2006.rds",sep="_")))
  return(non_con_model)
}


# where possible choose one of the modes for each nonconverged model -----------
get_inds <- function(i, nonconverged_rows){
  print(i)
  mod_obj <- get_model_obj(nonconverged_rows[i,], run_dir)
  # find indicies to discard
  disc_ind <- discard_chain_inds(mod_obj)
  return(disc_ind)
}


disc_inds <- map(1:(dim(nonconverged_rows)[1]), get_inds, nonconverged_rows)


# recalculate the assessment metrics for these models---------------------------
fixed_assess <- map(1:(dim(nonconverged_rows)[1]), function(i){
  print(i)
  disc_ind <- disc_inds[[i]]
  if (is.na(disc_ind)){
    return(NA)
  }
  mod_obj <- get_model_obj(nonconverged_rows[i,], run_dir)
  # get rates for non discarded chains
  fert_fit <- mod_obj$fert_fit
  stan_data <- mod_obj$stan_data
  rates_df <- get_rates(fert_fit, stan_data)
  rates_df %<>% mutate(num_Sim=as.numeric(gsub("V", "", Sim))) %>% 
    filter(!(num_Sim %in% disc_ind))
  disp <- as.matrix(fert_fit, "dispersion")[-disc_ind,]
  emp_rates_df <- get_empirical_rate_df("VV", nonconverged_rows[i,]$Country)
  nb_rates_df <- get_nb_rates(rates_df, disp, stan_data, emp_rates_df)
  
  # assess vs empirical
  forecast_jump_off <- (stan_data$country_config$last_year - 
                          stan_data$model_config$T_holdback)
  fore_score <- calc_score(nb_rates_df,
                           jump_off=forecast_jump_off,
                           axis=stan_data$model_config$time_dimension,
                           country=stan_data$country_config$country)
  fore_error <- calc_central_errors(nb_rates_df, forecast_jump_off, 
                                    stan_data$model_config$time_dimension,
                                    country=stan_data$country_config$country)
  
  
  # return new assessments from remaining chains
  return(cbind(nonconverged_rows[i,
                                 c("country", 
                                   "dispersion", 
                                   "param_type_1",
                                   "param_type_2")],
               fore_error$assess_sum_df, 
               Total_Score=fore_score$total_score,
               n_div=get_num_divergent(fert_fit)))
})


assess_out <- rbind(assessment_df %>% 
                      filter(n_div==0,
                             dispersion=="variable",
                             converged_rhat==T) %>%
                      select(country,
                             param_type_1,
                             param_type_2,
                             RMSE_rate, 
                             Q_90,
                             Q_50),
                    bind_rows(fixed_assess[!is.na(fixed_assess)]) %>%
                      select(country,
                             param_type_1,
                             param_type_2,
                             RMSE_rate, 
                             Q_90,
                             Q_50))

saveRDS(assess_out, "processed/assess_out.Rds")
