library(rstan)
library(tibble)
library(tidyr)
library(purrr)
library(magrittr)
library(yaml)
library(splines)
library(readr)
library(stringr)
library(dplyr)

source("R/data_functions.R")
source("R/stan_model_functions.R")
source("R/rate_reconstruction.R")

## Read configuations ----------------------------------------------------------
cmd_arg <- commandArgs(trailingOnly = T)
country <- cmd_arg[1]
config <- cmd_arg[2]

country_config <- read_yaml(file.path("config", "country",
                                      paste0(country, ".yaml")))

model_config <- read_yaml(file.path("config", "model",
                                    paste0(config, ".yaml")))


## Read in data ----------------------------------------------------------------
# Allow for either period or cohort models
lexis_prefix <- case_when(model_config$time_dimension=="Cohort"~"VV",
                          model_config$time_dimension=="Period" ~ "RR")

births_df <- readRDS(file.path("data", lexis_prefix, country, "births_hfd.rds"))
expos_df <- readRDS(file.path("data", lexis_prefix, country, "expos_hfd.rds"))


## make stan data --------------------------------------------------------------
stan_data <- make_stan_data(births_df, expos_df, model_config, country_config)

## construct model code (based on configuration) -------------------------------

fert_str <- get_fert_str(model_config$model,
                         model_config$dispersion, 
                         model_config$param_func)

inits <- function(){
  list(CEB= runif(stan_data$T_total, 0.5, 3),
       mode_age_sum=runif(stan_data$T_total, 50, 60),
       mode_age_gap=runif(stan_data$T_total, 5, 10),
       sd_age=matrix(runif(2 * stan_data$T_total,
                           4,8), 2, stan_data$T_total))
}

## fit the model ---------------------------------------------------------------
fert_model <- stan_model(model_code= fert_str, model_name = model_config$model)

fert_fit <- sampling(fert_model, data=stan_data, cores=3, chains=3, iter=2000,
                     init=inits,
                     control=list(max_treedepth=12,
                                  adapt_delta=0.95))

dir.create(model_config$out_dir,recursive = T)
saveRDS(fert_fit, file.path(model_config$out_dir, "test.rds"))


## assessment ------------------------------------------------------------------
rate <- get_rates(fert_fit, stan_data)

