library(tidyr)
library(ggplot2)
library(rstan)
library(dplyr)
library(magrittr)
library(scales)
library(tibble)
library(ggfan)
library(purrr)

source(file.path("R","data_functions.R"))
source(file.path("R","rate_reconstruction.R"))

results <- readRDS(file.path("mcmc_results",
                             "run_20190310_234223",
                             "GBRTENW",                       
                             paste0("mixture_gamma_density_",
                                    "weibull_density_variable_2006.rds")))
stan_data <- results$stan_data
fert_fit <- results$fert_fit
rm(results)
trash <-gc()

disp <- as.matrix(fert_fit, "dispersion")
rate_df <- get_rates(fert_fit, stan_data)
#nb_rate_df <- get_nb_rates(rate_df,disp, stan_data, emp_rates_df)

axis <- "Cohort"
rate_comp_df <- rate_df  %>% gather(Component, Rate, -!!sym(axis), -Age, -Sim) %>%
  mutate(Component=case_when(Component=="Comp_1" ~ "Component One",
                             Component=="Comp_2" ~ "Component Two",
                             Component=="Rate" ~ "Rate"))

get_quants <- function(df, probs=(1:99)/100){
  qs <- quantile(df$Rate, probs=probs)
  return(tibble(Quantile=probs, Rate=qs))
}

rate_comp_df  %<>% nest(-Age,-Cohort,-Component) %>% 
  mutate(Quant_df= map(data, get_quants)) %>% select(-data) %>% unnest()

saveRDS(rate_comp_df, "processed/GBRTENW_gamma_Weibull_rate_comp_df.rds")

