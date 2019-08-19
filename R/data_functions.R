
# Construct  a list of data to be passed into the stan fitting process,
# given birth and exposure data, and configuration for the model 
# and country specific configuratiosn
make_stan_data <- function(births_df, expos_df, model_config, country_config){
  
  births_df <- restrict_births_df(births_df, model_config, country_config)
  N <- dim(births_df)[1]
  
  # create data vector
  births <- births_df$Births %>% as.numeric()
  births <- round(births)
  
  # age index
  age_ind <- births_df$Age_i
  n_ages_obs <- length(unique(age_ind))
  
  # time index
  time_ind <- births_df$Time
  T_insample <- length(unique(births_df$Time))
  T_forecast <- model_config$T_forecast
  T_total <- T_insample + T_forecast
  
  # exposure
  expos <- get_exposure_mat(expos_df, model_config, country_config)
  expos[is.na(expos)] <- 0
  n_ages <- dim(expos)[1]
  
  # dispersion basis
  disps <- get_dispersion_basis(n_ages_obs, model_config$n_dispersion_basis - 2)
  
  ## dispersion index
  # This mostly maps to the age index
  # Except for the first and last ages
  births_df %<>% mutate(Disp_ind = ifelse(Age_i==1,
                                          n_ages_obs + Time,
                                          Age_i),
                        Disp_ind = ifelse(Disp_ind==n_ages_obs,
                                          n_ages_obs + T_insample + Time,
                                          Disp_ind))
  disp_ind <- births_df$Disp_ind
  
  # data gaps
  data_gap_low=country_config$first_data_age - model_config$minimum_age
  data_gap_high = model_config$maximum_age - country_config$last_data_age
  
  return(list(
    # lengths
    N=N,
    n_ages=n_ages,
    n_ages_obs=n_ages_obs,
    T_total=T_total,
    T_insample=T_insample,
    
    # data bounds
    minimum_age = model_config$minimum_age,
    data_gap_low = data_gap_low,
    data_gap_high = data_gap_high,
    
    # indexes
    age_ind=age_ind,
    time_ind=time_ind,
    disp_ind=disp_ind,
    
    # dispersion basis
    disp_basis=disps$disp_basis,
    n_basis=disps$n_basis,
    
    # data
    births=births,
    expos=expos,
    
    # config
    model_config=model_config,
    country_config=country_config
  ))
  
}

# Given raw HFD birth data, aggregate the extrapolated top and bottom ages.
restrict_births_df <- function(births_df, model_config, country_config){
  # extract a few index variables from the config files
  # this helps keep the code readable and the line length manageble!
  first_data_age <- country_config$first_data_age
  last_data_age <- country_config$last_data_age
  first_full_cohort <- country_config$first_full_cohort
  last_year <- country_config$last_year
  time_dim <- model_config$time_dimension
  first_full_cohort <- country_config$first_full_cohort
  
  # aggregating the ages below and above the observed data, 
  # where HFD extrapolates
  births_df %<>% mutate(Age=ifelse(Age >= first_data_age, 
                                   Age, first_data_age), 
                        Age=ifelse(Age <= last_data_age, 
                                   Age, last_data_age)) %>% 
    group_by(!!sym(time_dim), Age) %>% summarise(Births=sum(Births)) %>% 
    ungroup()
  
  ## Age index 
  births_df %<>% mutate(Age_i=Age - first_data_age + 1)
  
  # get rid of stuff falling outside our time horizons
  # create time index (either cohort or period)
  if(time_dim=="Cohort"){
    if(!("Cohort" %in% colnames(births_df))){
      stop("Cohort not present in data frame; conflict with specified 
           configuration")
    }
    births_df %<>% 
      mutate(Year=Cohort+Age) %>%
      filter(Cohort >= first_full_cohort,
             Year <= last_year - model_config$T_holdback) %>%
      mutate(Time=Cohort - min(Cohort) + 1) 
    
  } else if (run_config$time_dim=="period"){
    births_df %<>% 
      filter(Year >= max(data_config$first_full_year, 1950),
             Year <= last_year - model_config$T_holdback) %>%
      mutate(Time=Year - min(Year) + 1) 
  } else {
    stop("Time dimension must be cohort or period")
  }
  
  
  return(births_df)
}


# return the exposure matrix for model fitting
get_exposure_mat <- function(expos_df, model_config, country_config){
  last_year <- country_config$last_year
  time_dim <- model_config$time_dimension
  T_holdback <- model_config$T_holdback
  first_data_age <- country_config$first_data_age
  
  if(time_dim=="Cohort"){
    if(!("Cohort" %in% colnames(expos_df))){
      stop("Cohort not present in data frame; conflict with specified 
           configuration")
    }
    # Need three restrictions to avoid including those below first data age in
    # last year, e.g. 12 in 2016 -> 2004 cohort.
    expos <- expos_df %>% mutate(Year =Cohort + Age) %>%
      filter(Year <= country_config$last_year - T_holdback,
             Cohort <= country_config$last_year - T_holdback - first_data_age,
             Cohort >= country_config$first_full_cohort,
             Age >= model_config$minimum_age,
             Age <= model_config$maximum_age) %>%
      select(Age,Cohort,Exposure) %>% spread(Cohort, Exposure) %>%
      select(-Age) %>% as.matrix()
  } else if(time_dim=="Period"){
    expos <- expos_df %>%
      filter(Year <= country_config$last_year - T_holdback,
             Year >= max(1950, country_config$first_full_year)) %>%
      select(Age, Year, Exposure) %>% spread(Year, Exposure) %>%
      select(-Age) %>% as.matrix()
  } else {
    stop("Specified time dimension not understood. Must be 'Cohort' or 'Period")
  }
  return(expos)
}

# Create a spline basis function with which to fit dispersion
get_dispersion_basis <- function(n_age_obs, n_basis){
  # Make dispersion spline 
  disp_basis <- get_spline_basis(1:n_age_obs, 6)
  # force last two coefficients to be equal
  n_basis <- dim(disp_basis)[2]
  disp_basis[,n_basis-1] <- rowSums(disp_basis[,(n_basis-1):n_basis]) 
  disp_basis[,2] <- rowSums(disp_basis[,1:2])
  disp_basis<- disp_basis[,-c(1,n_basis)]
  n_basis <- n_basis - 2
  return(list(disp_basis=disp_basis, n_basis=n_basis))
}

# Create a basis function spline with equally spaced knots
get_spline_basis <- function(xx, n_interior_knots){
  knots_interior <- seq(min(xx), max(xx), length.out = n_interior_knots)
  knot_gap <- diff(knots_interior)[1]
  knots <-  c(min(xx) - (3:1* knot_gap), 
              knots_interior,
              max(xx) + 1:3*knot_gap)
  basis <- splineDesign(knots,xx)
  return(basis)
}

# Read in the empirical rates in the specifed country and format
get_empirical_rate_df <- function(hfd_code="VV", country="GBRTENW", prefix="."){
  births_df <- readRDS(file.path(prefix, "data", hfd_code, country, 
                                 "births_hfd.rds")) 
  expos_df <- readRDS(file.path(prefix, "data", hfd_code, country, 
                                "expos_hfd.rds")) 
  rate_df <- merge_data(expos_df, births_df) %>% as_tibble()
  
  if(!("Births" %in% colnames(rate_df))){
    rate_df %<>% mutate(Births=Total)
  }
  rate_df %<>% mutate(Rate= Births/Exposure)
  return(rate_df)
  
}

# put HFD birth and exposure data in single data frame and 
# calculate rates
merge_data <- function(expos_df, births_df){
  if (!("Births" %in% colnames(births_df))){
    births_df %<>% rename(Births=Total)
  }
  rate_df <- left_join(expos_df %>% select(-OpenInterval),
                       births_df %>% select(-OpenInterval)) %>% 
    mutate(Rate= Births/Exposure)
  return(rate_df)
}

