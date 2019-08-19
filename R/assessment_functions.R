
# For every sample of rates, simulate from the negative binomial distribution 
# with appropriate dispersion to obtain the correct posterior 
# predictive distribution for observed rates
get_nb_rates <- function(rate_df, disp, stan_data,emp_rates_df){
  first_data_age <- stan_data$country_config$first_data_age
  disp_df <- disp %>% as_tibble() %>% mutate(Sim=unique(rate_df$Sim))
  if(dim(disp)[2] > 1){
    disp_df %<>% gather(Age,Dispersion,-Sim) %>%
      mutate(Age=as.numeric(stringi::stri_extract(Age, regex="[0-9]+")) + 
                 first_data_age - 1)
  } else {
    colnames(disp_df)[1] <- "Dispersion"
  }
  
  rate_df %<>% left_join(disp_df) %>% 
    drop_na() %>%
    left_join(emp_rates_df %>% rename(Obs_rate=Rate)) %>%
    mutate(E_births = Rate*Exposure) %>%
    mutate(PP_births = rnbinom(n(), mu=E_births, size=Dispersion),
           PP_rate = PP_births/Exposure)
  return(rate_df) 
}

# calculate the pointwise log-score.
calc_score <- function(nb_rate_df, jump_off, axis,country=NULL){
  # filter for forecasts only
  assess_df <- nb_rate_df %>% filter(Year > jump_off,
                                     !is.na(Exposure),
                                     Exposure!=0) %>%
    mutate(Score = dnbinom(round(Births),
                           mu=E_births,
                           size = Dispersion,
                           log=T),
           Var_births = E_births + E_births**2 / Dispersion)
  
  if(axis=="Cohort"){
    col_name="Cohort"
  } else {
    col_name ="Year"
  }
  
  assess_sum_df <- assess_df %>% group_by(Age,!!sym(col_name)) %>% 
    summarise(Score=matrixStats::logSumExp(Score) - log(n()))
  
  total_score <- sum(assess_sum_df$Score)
  if(!is.null(country)){
    assess_df %<>% mutate(Country=country)
    assess_sum_df %<>% mutate(Country=country)
  }
  
  return(list(total_score=total_score, assess_df=assess_df,
              assess_sum_df=assess_sum_df, country=country))
}

# Calculate set of errors based on difference between 
# mean of predictive distribution and observation.
# Also calculate coverages
calc_central_errors <- function(assess_df, jump_off, axis, country=NULL){
  # log rates problematic (-inf), so ignored
  if(axis=="Cohort"){
    col_name="Cohort"
  } else {
    col_name ="Year"
  }
  assess_df %<>% filter(Year > jump_off,
                        !is.na(Exposure))
  # First, average over samples
  assess_df %<>% group_by(Age, !!sym(col_name)) %>% 
    summarise(Rate=mean(Rate),
              Log_rate=mean(log(Rate)),
              E_births=mean(E_births),
              Obs_rate=first(Obs_rate),
              Births=first(Births),
              Exposure=first(Exposure),
              Q = sum(PP_births <= Births)/n())
  
  # then, calculate the error measures pointwise
  assess_df %<>%
    mutate(Bias_rate = Rate-Obs_rate,
           AE_rate = abs(Bias_rate),
           ## calc PE of estimates rate, not obs, avoid / 0 !
           PE_rate = AE_rate/Rate, 
           SE_rate = (Bias_rate)**2,
           Bias_births=E_births-Births,
           AE_births=abs(Bias_births),
           PE_births=AE_births/E_births,
           SE_births=Bias_births**2)
  
  
  # average over each point
  assess_sum_df <- assess_df %>% ungroup() %>% 
    summarise(
      Q_90 = sum(as.numeric(Q>=0.05  & Q<=.95))/n(),
      Q_80 = sum(as.numeric(Q>=0.10  & Q<=.90))/n(),
      Q_50 = sum(as.numeric(Q>=0.25  & Q<=0.75))/n(),
      Bias_rate=mean(Bias_rate),
      AAE_rate = mean(AE_rate),
      APE_rate= mean(PE_rate),
      RMSE_rate = sqrt(mean(SE_rate)),
      Bias_births=mean(Bias_births),
      AAE_births=mean(AE_births),
      APE_births= mean(PE_births),
      RMSE_births=sqrt(mean(SE_births)))
  if(!is.null(country)){
    assess_df %<>% mutate(Country=country)
    assess_sum_df %<>% mutate(Country=country)
  }
  return(list(assess_df=assess_df, assess_sum_df=assess_sum_df))
}

# Has this set of samples converged, according the split rhat?
# Are there any divergent transition?
is_converged <- function(fert_fit){
  answer <- TRUE
  sum_fit <- summary(fert_fit)[[1]][1:get_num_parameters(fert_fit),]
  if(any(sum_fit[,10]>1.05)){
    print("rhat greater than 1.05")
    answer<-FALSE
  }
  
  n_div <- get_ndiv(fert_fit)
  if (n_div>0){
    print(paste0("Chains contains ",n_div, " divergent transitions"))
    answer=FALSE
  }
  
  return(answer)
}

# Return the number of divergent transitions
get_ndiv <- function(fert_fit){
  samp_pars <- get_sampler_params(fert_fit, inc_warmup=F)
  
  n_div <- map_dbl(samp_pars, function(mat){
    ind <- which(grepl("divergent",colnames(mat)))
    n_divergent <- sum(mat[,ind])
    return(n_divergent)
  })
  sum_n_div <- sum(n_div)
  return(sum_n_div)
}

# How many samples exceeded the maximum treedepth
get_n_max_tree_depths <- function(fert_fit){
  # from rstan::check_treedepth (which only outputs a message)
  # maxdepth from sampler args -> 10 by default
  max_depth <- fert_fit@stan_args[[1]]$control$max_treedepth
  if (is.null(max_depth)) {
    max_depth <- 10
  }
  treedepths <- rstan:::sampler_param_vector(fert_fit, "treedepth__")
  n <- sum(treedepths == max_depth)
  return(n)
}


get_summary_df <- function(res, country, axis, run_date){
  git_commit_time <- git2r::when(res$git_hash)
  git_commit_hash <- res$git_hash$sha
  git_commit_short_hash <- substr(res$git_hash$sha, 1, 7)
  elapsed_time <- max(rowSums(get_elapsed_time(res$fert_fit)))
  sum_fit <- summary(res$fert_fit)[[1]] %>% as_tibble() %>% drop_na()
  
  n_pars <- get_num_parameters(res$fert_fit)
  rhat <- sum_fit[1:n_pars, 10]
  n_eff <- sum_fit[1:n_pars, 9]
  n_div <-get_ndiv(res$fert_fit)
  
  energy_good <- utils::capture.output(rstan::check_energy(res$fert_fit), 
                                       type="message")
  energy_good <- grepl("no pathological", energy_good)
  
  # return 1 component if "single" is in model name
  n_components <- 1 + (1- grepl("single",
                                res$stan_data$model_config$model))
  mixed_components <- length(res$stan_data$model_config$param_fun) >1
  param_type_1 <- res$stan_data$model_config$param_fun[1]
  param_type_2 <- res$stan_data$model_config$param_fun[2]
  
  
  last_data_year <- res$stan_data$country_config$last_year
  forecast_jump_off <- last_data_year - res$stan_data$model_config$T_holdback
  if(axis=="Cohort"){
    start_time <- res$stan_data$country_config$first_full_cohort
  } else {
    start_time <- max(res$stan_data$country_config$first_full_year,1950)
  }
  
  #output_table
  out_df <- tibble(
    model_name=res$stan_data$model_config$model,
    country=country,
    axis=axis,
    dispersion= res$stan_data$model_config$dispersion,
    n_components = n_components,
    mixed_components = mixed_components,
    param_type_1 = param_type_1,
    param_type_2 = param_type_2,
    converged_rhat = max(rhat)<1.05,
    n_large_rhat = sum(rhat>1.05),
    min_neff = min(n_eff),
    runtime_min = elapsed_time/60,
    energy_good=energy_good,
    treedepth_good= get_n_max_tree_depths(res$fert_fit)==0,
    n_max_treedepth= get_n_max_tree_depths(res$fert_fit),
    n_div=n_div,
    git_commit_time=git_commit_time,
    git_commit_hash=git_commit_hash,
    git_commit_short_hash=git_commit_short_hash,
    run_date=run_date,
    n_holdback=res$stan_data$model_config$T_holdback,
    forecast_jump_off=forecast_jump_off,
    last_data_year=res$stan_data$country_config$last_year
  ) 
  return(out_df)
}

# Produce a summary dataframe of the posterior sample
get_sum_fit_df <- function(fert_fit , chain="all"){
  npars <- get_num_parameters(fert_fit)
  if (chain=="all"){
    sum_fit <- summary(fert_fit)[[1]][1:npars,]
  } else {
    sum_fit <- summary(fert_fit)[[2]][1:npars,,chain]
  }
  
  pars <- rownames(sum_fit)[1:npars]
  sum_fit_df <- sum_fit %>% 
    as_tibble() %>%
    mutate(Parameter=pars,
           Family=gsub("\\[[0-9]{1,},*[0-9]*\\]",
                       "",
                       Parameter, perl=T))
  return(sum_fit_df)
  
}

# How many parameters are directly sampled in the model?
get_num_parameters <- function(fert_fit){
  split_adapt_strings <- stringr::str_split(rstan::get_adaptation_info(fert_fit), "#")
  n_pars <- stringr::str_split(split_adapt_strings[[1]][5], ", ") %>% 
    unlist() %>%
    length()
  return(n_pars)
}

# What were the stepsizes
get_stepsizes <- function(fert_fit){
  split_adapt_strings <- stringr::str_split(rstan::get_adaptation_info(fert_fit), "#")
  stepsize <- map_dbl(split_adapt_strings, 
                      function(x) as.numeric(stringr::str_extract(x[3],"[0-9]+\\.[0-9]+")))
  return(stepsize)
}

# Read rstan diagnostic files
read_diag <- function(i, file_name, diag_dir="results"){
  diag_df <- read.csv(file.path(diag_dir,
                                paste0(file_name, "_", i, ".csv")),
                      comment.char="#") %>% as_tibble() %>%
    mutate(iter=1:n(), Chain=i)
  return(diag_df)
}

# Calculate the loo diagnostic
get_loo <- function(fert_fit, stan_data, do_ind=F){
  if(do_ind){
    inds <- stan_data$age_ind > 3 & stan_data$age_ind <33
  } else {
    inds <- rep(TRUE, length(stan_data$age_ind))
  }
  
  ll <- extract_log_lik(fert_fit, merge_chains = F)[,,inds]
  r_eff <- relative_eff(ll)
  loo_obj <- loo(ll, r_eff=r_eff)
  return(loo_obj)
}

# return a dataframe with pointwise pareto K values, with appropriate indexes
loo_diagnostics <-function(loo_obj, stan_data, inds=NULL){
  
  if (is.null(inds)){
    inds <- rep(TRUE, stan_data$N)
  }
  first_age <- stan_data$country_config$first_data_age
  first_cohort <- stan_data$country_config$first_full_cohort
  loo_df <- tibble(Pareto_k = loo_obj$diagnostics$pareto_k[inds],
                   Age=stan_data$age_ind[inds] + first_age - 1,
                   Cohort = stan_data$time_ind[inds] + first_cohort - 1)
  return(loo_df)                   
}


