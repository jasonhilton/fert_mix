
# Given the fitted fertility model (of class stanfit), and the data used to fit 
# the model
# extract the mixture time series from the fitted stan object as a data frame
get_mixture <- function(fert_fit, stan_data){
  if( stan_data$model_config$time_dimension== "Cohort"){
    axis <- "Cohort"
    start_time <- stan_data$country_config$first_full_cohort
  } else {
    axis <- "Year"
    start_time <- max(1950, stan_data$country_config$first_full_cohort)
  }
  
  logit_mixture <- get_time_series(fert_fit, stan_data, param_name="logit_mixture",
                   col_name = "Logit_mixture",initial_value = start_time,
                   time_dim = axis)
  mixture <- logit_mixture %>% mutate(Mixture = plogis(Logit_mixture)) %>% 
    select(-Logit_mixture)

  return(mixture)
}


# Given the fitted fertility model (of class stanfit), and the data used to fit 
# the model
# extract the time series named "param_name". 
# In the return dataframe, the columns relating to the parameter value and the 
# time index will be named with the value of the string `col_name` and 
# `time_dim` respectively, with the latter index starting from `initial_value`
get_time_series <- function(fert_fit, stan_data, param_name,
                            col_name, initial_value,
                            time_dim){
  # if a name for the value column in returned dataframe is specified, use this
  # Otherwise, use the (stan) param_name.
  if( is.null(col_name)){
    col_name <- param_name
  }
  var_matrix <- as.matrix(fert_fit, param_name)
  
  long_df <- var_matrix %>% t() %>% tibble::as_tibble() %>% 
    dplyr::mutate(!!time_dim:= initial_value +  1:stan_data$T_total - 1) %>%
    tidyr::gather(Sim, !!col_name, -!!time_dim)
  
  return(long_df)
}


# Returns a function corresponding to a parametric density with name `func_name`
# Must be one of "weibull_density", "gamma_density" or "hadwiger density
funct_grabber <- function(func_name){
  if (func_name=="weibull_density"){
    out_func <- function(x, mu, sd){
      upper_quartile <- mu + sd
      k <- log(2) / (log(upper_quartile) - log(mu))
      lam = mu / ((log(2))**(1.0/k))
      return(dweibull(x = x, shape = k, scale=lam))
    }
  } else if(func_name=="gamma_density"){
    out_func <- function(x, mu, sd){
      rate = (mu + sqrt(mu^2 + 4*sd^2)) / (2*sd^2)
      shape = mu * rate + 1
      return(dgamma(x = x, rate = rate, shape = shape))
    }
  } else if (func_name=="hadwiger_density"){
    out_func <- function(x, location, scale){
      log_H = log(location) - 0.5*(log(2)) - log(scale);
      part1 = log_H - log(location) - 0.5 * log(pi);
      part2 = 1.5  * log(location) - 1.5 * log(x);
      part3 = - (exp(log_H)^2)  * (location/x + x/location - 2);
      
      return (exp(part1 + part2 + part3));
    }
    
  }
  return(out_func)
}

# Given the fitted fertility model (of class stanfit), and the data used to fit 
# the model
# return a data frame containing the location parameters of the parametric 
# components.
get_modes <- function(fert_fit, stan_data){
  if( stan_data$model_config$time_dimension== "Cohort"){
    axis <- "Cohort"
    start_time <- stan_data$country_config$first_full_cohort
  } else {
    axis <- "Year"
    start_time <- max(1950, stan_data$country_config$first_full_cohort)
  }
  mode_age_gap <- get_time_series(fert_fit, stan_data, 
                                  param_name="mode_age_gap",
                                  col_name = "Mode_gap",
                                  initial_value = start_time,
                                  time_dim = axis)
  mode_age_sum <- get_time_series(fert_fit, stan_data, 
                                  param_name="mode_age_sum",
                                  col_name = "Mode_sum",
                                  initial_value = start_time,
                                  time_dim = axis)
  
  mode_df <- left_join(mode_age_gap, mode_age_sum)
  mode_df %<>% mutate(Mode_1 = 0.5*(Mode_sum - Mode_gap),
                      Mode_2 = 0.5*(Mode_sum + Mode_gap)) %>%
    select(-Mode_sum, -Mode_gap)
  return(mode_df)
}

# Given the fitted fertility model (of class stanfit), and the data used to fit 
# the model
# return a data frame containing the spread parameters of the parametric 
# components.
get_sd <- function(fert_fit, stan_data){
  if( stan_data$model_config$time_dimension== "Cohort"){
    axis <- "Cohort"
    start_time <- stan_data$country_config$first_full_cohort
  } else {
    axis <- "Year"
    start_time <- max(1950, stan_data$country_config$first_full_cohort)
  }
  sd_df <- get_2d_stan_df(fert_fit,"sd_age","Component", axis,
                          conversion_2_func = function(x) x+start_time -1,
                          value_name = "Sd")
  
  sd_df <- spread(sd_df, Component, Sd) %>% rename(Sd_1 = `1`,
                                                   Sd_2 = `2`)
  return(sd_df)
}

# Given the data used to fit the ferility model and the various timeseries
# return a data frame containing the two components making up the fertility 
# schedule.
get_raw_components <- function(stan_data, mode_df, sd_df, func_1,
                               func_2){
  # func_1 is LOWEST function here
  if( stan_data$model_config$time_dimension== "Cohort"){
    axis <- "Cohort"
    
  } else {
    axis <- "Year"
  }
  
  comp_df <- left_join(mode_df, sd_df)
  age_min <- stan_data$minimum_age
  n_ages <- stan_data$n_ages
  comp_df %<>% mutate(Age=list(1:n_ages)) %>% unnest()
  comp_df %<>% mutate(Comp_1 = func_1(Age - 0.5, 
                                      Mode_1 - age_min, Sd_1),
                      Comp_2 = func_2(Age - 0.5, 
                                      Mode_2 - age_min, Sd_2))
  comp_df %<>% mutate(Age=Age + age_min - 1)
  return(comp_df %>% select(!!axis, Age, Sim, Comp_1, Comp_2))
  
}

# Given the fitted fertility model (of class stanfit), and the data used to fit 
# the model
# Return a dataframe of posterior rate estimates from the fitted model
get_rates <- function(fert_fit,stan_data){
  CEB_df <- get_CEB(fert_fit,stan_data)
  mode_df <- get_modes(fert_fit, stan_data)
  sd_df <- get_sd(fert_fit, stan_data)
  mixture_df <- get_mixture(fert_fit, stan_data)
  rate_df <- get_comps(stan_data, CEB_df, mode_df ,sd_df, mixture_df)
  return(rate_df)
}

# return components making up the fertility curve, given the timeseries of the 
# model paraemeters
get_comps <- function(stan_data, CEB_df, mode_df, sd_df, mixture_df){
  if( stan_data$model_config$time_dim== "Cohort"){
    axis <- "Cohort"
  } else {
    axis <- "Year"
  }
  
  func_low <- funct_grabber(stan_data$model_config$param_func[1])
  func_high <- funct_grabber(stan_data$model_config$param_func[2])
  
  raw_comps <- get_raw_components(stan_data, 
                                  mode_df, 
                                  sd_df,
                                  func_low, 
                                  func_high)
  
  comp_df <- left_join(raw_comps, CEB_df) %>% left_join(mixture_df) %>%
    mutate(Comp_1 = CEB*(Mixture)*Comp_1,
           Comp_2 = CEB*(1-Mixture)*Comp_2,
           Rate = Comp_1 + Comp_2) %>% select(Age,!!axis, Sim, 
                                              Comp_1, Comp_2, Rate)
  return(comp_df)
}

# Given the fitted fertility model (of class stanfit), and the data used to fit 
# the model
# return a data frame containing the fertility summary measure (Children ever born)
get_CEB <- function(fert_fit, stan_data){
  if( stan_data$model_config$time_dim== "Cohort"){
    axis <- "Cohort"
    sum_name <- "CEB"
    start_time <- stan_data$country_config$first_full_cohort
  } else {
    axis <- "Year"
    sum_name <- "TFR"
    start_time <- max(1950, stan_data$country_config$first_year)
  }
  
  if ("CEB" %in% fert_fit@model_pars){
    CEB_df <- get_time_series(fert_fit, stan_data, param_name="CEB",
                              col_name = sum_name,
                              initial_value = start_time,
                              time_dim = axis)  
  } else {
    stop("Can't find CEB in model")
  }
  
  return(CEB_df)
}

# extract a generic matrix from a stan model into a tidy data frame,
# potential with transformed indicies.
get_2d_stan_df  <- function(model_fit, param_name,
                            index_1_name= "Index_1", 
                            index_2_name= "Index_2",
                            conversion_1_func=identity,
                            conversion_2_func=identity,
                            value_name ="Value"){
  if(class(model_fit)!="stanfit"){
    stop("model_fit must be of class `stanfit`")
  }
  par_mat <- as.matrix(model_fit, param_name)
  
  par_names <- colnames(par_mat)
  
  if(max(stringr::str_count(par_names,"\\,")) > 1){
    stop("Parameter has more than 2 dimensions")
  }
  if(min(stringr::str_count(par_names, "\\,")) <1 ){
    stop("Parameter has fewer than 2 dimensions")
  }
  ind_1 <- stringr::str_extract_all(par_names,
                                    # match all numbers followed by a comma
                                    "[0-9]+(?=\\,)",
                                    simplify = T) %>% 
    as.numeric()
  ind_2 <- stringr::str_extract_all(par_names,
                                    # match all numbers preceded by a comma
                                    "(?<=\\,)[0-9]+",
                                    simplify = T) %>% 
    as.numeric()
  
  par_df <- par_mat %>% t() %>% tibble::as_tibble() %>%
    mutate(!!index_1_name:=conversion_1_func(ind_1),
           !!index_2_name:=conversion_2_func(ind_2)) %>%
    gather(Sim, !!value_name, -!!index_1_name, -!!index_2_name)
  return(par_df)
}

