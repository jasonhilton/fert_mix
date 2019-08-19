# ------------------------------------------------------------------------------

# list of strings to be inserted into stan template files, specifying different
# dispersion configurations.

dispersions <- list(
  constant = list(
    dispersion_data = "",
    dispersion_parameters = "real log_dispersion;",
    dispersion_transform = paste("real dispersion;",
                                 "dispersion=exp(log_dispersion);", sep="\n"),
    dispersion_sum_dec = "",
    dispersion_sums = "",
    lik_dist="neg_binomial_2_lpmf",
    dispersion_lik=", dispersion",
    dispersion_lik_ind="",
    dispersion_distribution= ""
  ),
  variable = list(
    dispersion_data = paste("int n_basis;",
                            "int disp_ind[N];",
                            "matrix[n_ages_obs, n_basis] disp_basis;",
                            sep="\n"),
    dispersion_parameters = paste("vector[n_basis] beta_disp;",
                                  "real<lower=0> sigma_disp;",
                                  sep="\n"),
    dispersion_transform = paste("vector[n_ages_obs] dispersion;",
                                 "dispersion = 1 + exp(disp_basis * beta_disp);",
                                 sep="\n"),
    dispersion_sum_dec = paste(
      "vector[n_ages_obs + 2 * T_insample] dispersion_star;",
      "dispersion_star[1:n_ages_obs] = dispersion;",
      sep="\n"),
    dispersion_sums = "for (t in 1:T_insample){ \n
    // Sum up the bottom \n
    dispersion_star[n_ages_obs + t] = (dispersion[1] *  \n
                    sum(E_births[1:(data_gap_low + 1),t]) ^ 2 / \n
                    dot_self(E_births[1:(data_gap_low + 1),t])); \n
    // Sum up the top \n
    dispersion_star[n_ages_obs + T_insample + t] = (dispersion[n_ages_obs] * \n
                    sum(E_births[(n_ages_obs + data_gap_low):n_ages,t]) ^ 2 / \n
                    dot_self(E_births[(n_ages_obs + data_gap_low):n_ages,t])); \n
  }\n",
    lik_dist="neg_binomial_2_lpmf",
    dispersion_lik=", dispersion_star",
    dispersion_lik_ind="[disp_ind[i]]",
    dispersion_distribution= paste(
      "beta_disp[1] ~ normal(0,5);",
      "beta_disp[2:n_basis] ~ normal(beta_disp[1:(n_basis-1)],sigma_disp*0.1);",
      "sigma_disp ~ normal(0,5);",sep="\n")
  ),
  poisson = list(
    dispersion_data = "",
    dispersion_parameters = "",
    dispersion_transform = "",
    dispersion_sum_dec = "",
    dispersion_sums = "",
    lik_dist="poisson_lpmf",
    dispersion_lik="",
    dispersion_lik_ind="",
    dispersion_distribution= ""
  ),
  var_unsmooth = list(
    dispersion_data = "int disp_ind[N];",
    dispersion_parameters = paste("vector<lower=0>[n_ages_obs] log_dispersion;",
                                  "real<lower=0> sigma_disp;", sep="\n"),
    dispersion_transform = paste("vector[n_ages_obs] dispersion;",
                                 "dispersion = 1.0 + exp(log_dispersion);",
                                 sep="\n"),
    dispersion_sum_dec = paste(
      "vector[n_ages_obs + 2 * T_insample] dispersion_star;",
      "dispersion_star[1:n_ages_obs] = dispersion;",
      sep="\n"),
    dispersion_sums = "for (t in 1:T_insample){ \n
    // Sum up the bottom \n
    dispersion_star[n_ages_obs + t] = (dispersion[1] *  \n
                    sum(E_births[1:(data_gap_low + 1),t]) ^ 2 / \n
                    dot_self(E_births[1:(data_gap_low + 1),t])); \n
    // Sum up the top \n
    dispersion_star[n_ages_obs + T_insample + t] = (dispersion[n_ages_obs] * \n
                    sum(E_births[(n_ages_obs + data_gap_low):n_ages,t]) ^ 2 / \n
                    dot_self(E_births[(n_ages_obs + data_gap_low):n_ages,t])); \n
  }\n",
    lik_dist="neg_binomial_2_lpmf",
    dispersion_lik=", dispersion_star",
    dispersion_lik_ind="[disp_ind[i]]",
    dispersion_distribution= paste(
      "log_dispersion[2:n_ages_obs] ~ normal(log_dispersion[1:(n_ages_obs-1)], sigma_disp);",
      "log_dispersion[1] ~ normal(0,4);",
      "sigma_disp ~ normal(0,1);",sep="\n")
  )
)


## Functions to create the stan model file--------------------------------------
  
# given the name of a stan model file, a specification of overdispersion and 
# a list of maximum two parametric functions
# read in the template and construct the full stan model
get_fert_str <- function(model_name, dispersion, param_funcs){
  fert_template <- read_file(file.path("stan", paste0(model_name, ".stan")))
  insert_list <- dispersions[[dispersion]]
  if (length(param_funcs)==2){
    param_fun <- read_file(file.path("stan",
                                     paste0("param_funcs", ".stan")))
    #                                 paste0("all_funcs_diff", ".stan")))
    insert_list <- c(insert_list, list(fun_name_1=param_funcs[1],
                                       fun_name_2=param_funcs[2]))
  } else {
    param_fun <- read_file(file.path("stan", "param_funs",
                                     paste0(param_funcs, ".stan")))  
  }
  insert_list$param_fun <- param_fun
  fert_str <- str_interp(fert_template, insert_list)
  return(fert_str)
}


# Create a data frame containing a set of model specifications (one per row)
# base settings are take from a config file, while all combinations of specified
# dispersions, model names, parametric functions and countries are constructed
# Optionally, functions with the same parametric function in the first and second
# position can be excluded
gen_runs_df <- function(base_config_file, param_types,
                        model_names, countries,
                        dispersion, exclude_dupls=F){
  par_grid <- expand.grid(param_type_1=param_types,
                          param_type_2=param_types,
                          model_name=model_names,
                          country=countries,
                          dispersion=dispersion,
                          stringsAsFactors = F)
  if (exclude_dupls){
    par_grid <- par_grid[par_grid$param_type_1 != par_grid$param_type_2,]  
  }
  
  par_grid$base_config <- base_config_file
  return(par_grid)
}

# given a base configuration (a nmaed list), and a model specification,
# create another named list, consisting of the base configuration with the values
#  of elements also appearing in run specification replacing existing values.
replace_base_config <- function(run_spec, base_config){
  config <- base_config
  for(par in names(run_spec)){
    config[par] = run_spec[par]
  }
  config$param_func <- c(config$param_type_1, config$param_type_2)
  return(config)
}