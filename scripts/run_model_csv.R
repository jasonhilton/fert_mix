library(rstan)
library(tibble)
library(tidyr)
library(purrr)
library(magrittr)
library(yaml)
library(splines)
library(readr)
library(stringr)
library(git2r)
library(loo)
library(dplyr)

source("R/data_functions.R")
source("R/stan_model_functions.R")
source("R/rate_reconstruction.R")
source("R/assessment_functions.R")

## Read csv --------------------------------------------------------------------

cmd_arg <- commandArgs(trailingOnly = T)
csv_line <- cmd_arg[1]
csv_file <- cmd_arg[2]

if (is.na(csv_file)){
  csv_file <- tail(list.files("run_specs"),1)
} 

csv_name <- unlist(strsplit(csv_file, "\\."))[1]

run_specs <- read.csv(file.path("run_specs", csv_file),
                      stringsAsFactors = F)

run_spec <- run_specs[csv_line,]
config <- unlist(run_spec["base_config"])

country <- run_spec$country 


country_config <- read_yaml(file.path("config", "country",
                                      paste0(country, ".yaml")))

model_config <- read_yaml(file.path("config", "model",
                                    paste0(config, ".yaml")))

model_config <- replace_base_config(run_spec, model_config)

cat(paste0("** config file name = ", config, "\n"))
cat(paste0("** csv file name = ", csv_file, "\n"))
cat(paste0("** csv run = ", csv_line, "\n"))
cat(paste0("** model_name = ", model_config$model, "\n"))
cat(paste0("** dispersion = ", model_config$dispersion, "\n"))
cat(paste0("** param_types = ", model_config$param_func, "\n"))
cat(paste0("** country = ", country, "\n"))
cat(paste0("** time ="), format(Sys.time(), "%Y%m%d_%H%M%S"), "\n")


model_name <- model_config$model
param_type <- model_config$param_func
forecast_jump_off <- country_config$last_year - model_config$T_holdback
dispersion <- model_config$dispersion


## Read in data ----------------------------------------------------------------
# Allow for either period or cohort models
lexis_prefix <- case_when(model_config$time_dimension=="Cohort"~"VV",
                          model_config$time_dimension=="Period" ~ "RR")

births_df <- readRDS(file.path("data", lexis_prefix, country, "births_hfd.rds"))
expos_df <- readRDS(file.path("data", lexis_prefix, country, "expos_hfd.rds"))


## make stan data --------------------------------------------------------------
stan_data <- make_stan_data(births_df, expos_df, model_config, country_config)

## construct model code (based on configuration) -------------------------------
fert_str <- get_fert_str(model_config$model, model_config$dispersion, 
                         model_config$param_func)

# set some reasonable starting values
inits <- function(){
  list(CEB= runif(stan_data$T_total, 0.5, 3),
       mode_age_sum=runif(stan_data$T_total, 50, 60),
       mode_age_gap=runif(stan_data$T_total, 5, 10),
       sd_age=matrix(runif(2 * stan_data$T_total,
                           4, 8), 2, stan_data$T_total),
       logit_mixture=runif(stan_data$T_total, 0.3,0.7))
}

## fit the model ---------------------------------------------------------------
fert_model <- stan_model(model_code= fert_str, model_name = model_config$model)

fert_fit <- sampling(fert_model, data=stan_data, cores=4, chains=4, iter=10000,
                     thin=4,
                     init=inits,
                     control=list(max_treedepth=12, adapt_delta=0.95))

dir.create(file.path("mcmc_results", csv_name, 
                     country), recursive=T)

result <-list(stan_data=stan_data,
              fert_fit=fert_fit,
              git_hash=commits(repository("."))[[1]])


cat("saving result ... \n")
cat(paste0("** time = "), format(Sys.time(), "%Y%m%d_%H%M%S"), "\n")
res_file_name <- paste0(model_name, "_",
                        paste(param_type, collapse="_"), "_",
                        dispersion, "_",
                        forecast_jump_off, ".rds")
saveRDS(result,
        file.path("mcmc_results", csv_name, country,
                  res_file_name))


## assessment ------------------------------------------------------------------
rate <- get_rates(fert_fit, stan_data)
disp <- as.matrix(fert_fit, "dispersion")
emp_rates_df <- get_empirical_rate_df(lexis_prefix,country)
nb_rates_df <- get_nb_rates(rate, disp, stan_data, emp_rates_df)


fore_score <- calc_score(nb_rates_df,
                         jump_off=forecast_jump_off,
                         axis=model_config$time_dimension,
                         country=country)
fore_error <- calc_central_errors(nb_rates_df, forecast_jump_off, 
                                  model_config$time_dimension,
                                  country=country)

run_date <- gsub("run_", "", csv_name)
sum_df <- get_summary_df(result, country, model_config$time_dimension, run_date)
sum_df <- cbind(sum_df, fore_error$assess_sum_df)
sum_df$total_score <- fore_score$total_score

ll <- extract_log_lik(fert_fit,merge_chains = F)
# loo_obj <- loo(ll[,chain_inds,],
#                r_eff = relative_eff(ll[,chain_inds,]), cores=4)
loo_obj <- loo(ll,
               r_eff = relative_eff(ll), cores=4)
sum_df$looic <- loo_obj$estimates[3,1]
sum_df$looic_se <- loo_obj$estimates[3,2]
sum_df$n_high_pareto_k <- sum(loo_obj$diagnostics$pareto_k>0.7)

## save assessements -----------------------------------------------------------
out_path <- file.path("mcmc_results", csv_name)
saveRDS(sum_df,
        file.path(out_path, 
                  paste0("assessment", "_",
                         csv_line, ".rds")))


## copy files ------------------------------------------------------------------
if (file.exists("backup_copy_file.sh")){
  # copy assessment file
  command <- paste0("backup_copy_file.sh ", out_path,
                    " assessment", "_", csv_line, ".rds")
  system(command)
}

if (file.exists("backup_copy_file.sh")){
  # copy result_file
  command <- paste0("backup_copy_file.sh ", out_path,"/", country, 
                    " ", res_file_name)
  system(command)
}

cat("teleiwsa bre! \n")
cat(paste0("** time = "), format(Sys.time(), "%Y%m%d_%H%M%S") ,"\n")

