source("R/data_functions.R")
source("R/stan_model_functions.R")

cmd_args <- commandArgs(trailingOnly = T)

config_file <- cmd_args[1]
if( is.na(config_file)){
  config_file <- "base_hb"
}

dir.create("run_specs")

countries <- c("GBRTENW", "FRATNP", "SWE", "USA")

dispersion <- c("variable")
model_names <- c("mixture")
param_types <- c("gamma_density", "weibull_density", "hadwiger_density")

runs_df <- gen_runs_df(config_file, param_types, model_names, countries,
                       dispersion)

time_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

readr::write_csv(runs_df, 
                 path=file.path("run_specs",
                                paste0("run_",
                                       time_stamp,
                                       ".csv")))
