

# create loo object given a specific file
#
#
create_loo <- function(mod_file, country, run_dir){
  file_path <- file.path(run_dir, country, mod_file)
  mod_obj <- readRDS(file_path)
  fert_fit <- mod_obj$fert_fit
  ll <- extract_log_lik(fert_fit, merge_chains = F)
  if (get_ndiv(fert_fit)>0){
    return(NA)
  }
  rhats <- summary(mod_obj$fert_fit)[[1]][,10]
  if (sum(rhats>1.05)>0){
    disc <- discard_chain_inds(mod_obj)
    n_iter <- dim(ll)[1] * dim(ll)[2]
    neg_chains <- chains_from_inds(disc,n_iter)
    if(is.na(neg_chains)){
      return(NA)
    }
    ll <- ll[,-neg_chains,]
  }
  
  r_eff <- relative_eff(ll)
  loo_obj <- loo(ll, r_eff=r_eff)
  return(loo_obj)
}

# calculate and save out loos for each country
compute_country_loos <- function(country,run_dir){
  mod_files <- list.files(file.path(run_dir, country))
  mod_files <- mod_files[stringr::str_detect(mod_files, "variable")]
  message("Computing loos...")
  walk(mod_files, function(mod_file, country, run_dir){
    message("stanfit file = ",mod_file)
    param_types <- str_extract_all(mod_file, "[a-z]+(?=_density)")
    param_types <- paste(unlist(param_types), collapse="_")
    loo_obj <- create_loo(mod_file, country, run_dir)
    dir.create(file.path("processed", country),recursive = T)
    saveRDS(loo_obj, file.path("processed", country,
                               paste0(param_types, "_loo.rds")))
  },
  country,
  run_dir)
}


##  given a list of data frame with a list-column of ll indexes to be 
# dropped from the model fitting and corresponding
# column splitting these indexes into <=3 k-fold groups,
# This function returns the a df with the ll indexes only from the `ind`-th group.
get_group_df <- function(ind, drop_df){
  group_ind <- drop_df %>%
    mutate(inds = map2_chr(group_indexes,drop_indexes,
                           function(group, drop){
                             paste(drop[group==ind], collapse="_")
                           }))
  return(group_ind)
}

# given a set of indicies to drop, and the corresponding 
# set of covarianes
# work out a non-terrible set of cv-fold/groups
# in order to estimate loo. 
generate_groups <- function(drop_inds, cohort_inds){
  if (length(drop_inds)==0 || is.na(drop_inds)){
    return(NA)
  }
  dist_mat <- cohort_inds %>% filter(ind %in% drop_inds) %>% 
    select(cohort,ages) %>% dist()
  n_drops <- length(drop_inds)
  dist_mat <- as.matrix(dist_mat)
  # don't consider distances from yourself as a minimum
  diag(dist_mat) <- 1000
  if (n_drops==1){
    return(1)
  } else if(n_drops==2) {
    return(c(1,2))
  } else if (n_drops==3){
    return(c(1,2,3))
  } else {
    # test k-fold cross validation groups to maximise the smallest within between
    # groups
    # 3 groups
    # 20 trials
    max_dist <- 0
    rand_split <- rep(c(1,2,3), ceiling(n_drops/2))[1:n_drops]
    best_split <- rand_split
    for (j in 1:100){
      # the smallest distance between any member of different groups
      this_try_min <- min(dist_mat[rand_split==1,rand_split==1])  
      this_try_min <- min(c(this_try_min,
                            as.numeric(dist_mat[rand_split==2,rand_split==2])))
      this_try_min <- min(c(this_try_min,
                            as.numeric(dist_mat[rand_split==3,rand_split==3])))
      if (this_try_min>max_dist){
        max_dist <- this_try_min
        best_split <- rand_split
      }
      rand_split <- kfold_split_random(3,n_drops)
    }
    if(max_dist <2){
      warning("Minimum distance within cv-fold is small!")
    } else{
      print (paste0("Minimum euclidean distance within cd-folds",
                    " in age-cohort space is ", max_dist))
    }
    return(best_split)
  }
}



# return a data frame for a specified country where each row corresponds to a 
# specification for one cv run to be carried out.
# The cv runs are chosen based on leaving out observations which have high 
# pareto-k diagnostics in the approximate PSIS loo calcutaions.
get_country_lko_run_list <- function(country,run_dir, pattern="variable"){
  
  # get information about each run.
  mod_files <- list.files(file.path(run_dir, country), pattern=pattern)
  model_index <- tibble(Model=paste0("model",1:9),file_name=mod_files) %>%
    mutate(param_type_1=map_chr(file_name, 
                                function(x) str_split(x, "_", 
                                                      simplify=T)[2])) %>% 
    mutate(param_type_2=map_chr(file_name, 
                                function(x) str_split(x, "_", 
                                                      simplify=T)[4])) %>%
    mutate(country=country)
  
  # read in loos
  loos <- map(seq(length(mod_files)), function(i){
    file_name <- paste(model_index$param_type_1[i], model_index$param_type_2[i],
                       "loo.rds", 
                       sep="_")
    return(readRDS(file.path("processed", country, file_name)))
  })
  
  # create list of indices to be rerun
  drop_ll <- map(loos, 
                 function(loo_obj){
                   if(is.na(loo_obj)){
                     return(NA)
                   } else {
                     which(loo_obj$diagnostics$pareto_k > 0.7)
                   }
                 }
  )
  
  # We want to create CV groups, where each group contains observations that 
  # are roughly independent. This is done by maximising the minimum l2 distance
  # between groups through random trials.
  # covariate inds
  
  stan_data <- readRDS(file.path(run_dir, country, mod_files[[1]]))$stan_data
  variate_inds <- tibble(cohort=stan_data$time_ind, 
                         ages = stan_data$age_ind,
                         ind=1:stan_data$N)
  
  group_indexes <- map(drop_ll,generate_groups, variate_inds)
  
  drop_df <- model_index %>% mutate(group_indexes=group_indexes,
                                    drop_indexes=drop_ll)
  
  
  new_runs <- map_df(1:3, get_group_df, drop_df) %>% 
    select(-group_indexes, -drop_indexes) %>% filter(inds!="" & inds!="NA") %>%
    mutate(dispersion=pattern, country=country)
  
  return(new_runs)
}

get_correct_loo <- function(run_row, run_spec_lko, run_dir){
  country <- run_row$country
  
  cntry_run_dir <- file.path(run_dir, country)
  pt1 <- gsub("_density", "", run_row$param_type_1)
  pt2 <- gsub("_density", "", run_row$param_type_2)
  message("Model = ", pt1, " / ", pt2)
  loo_obj <- readRDS(file.path("processed",country, 
                               paste(pt1,
                                     pt2,
                                     "loo.rds", sep="_")))
  if(is.na(loo_obj)){
    return(tibble(country=cntry, 
                  param_type_1=pt1,
                  param_type_2=pt2,
                  looic=NA,
                  looic_se=NA))
  }
  
  n_high <- sum(loo_obj$diagnostics$pareto_k>0.7)
  if (n_high==0){
    return(tibble(country=cntry, 
                  param_type_1=pt1,
                  param_type_2=pt2,
                  looic=loo_obj$estimates[3,1],
                  looic_se=loo_obj$estimates[3,2]))
  } else {
    cnrty_run_dir <- file.path(run_dir,country)
    
    out_df <- run_spec_lko %>% filter(country==cntry,
                                      param_type_1==pt1,
                                      param_type_2==pt2) %>% 
      mutate(lpd_df = map2(Index, inds, get_lpd, cntry_run_dir)) %>% unnest()
    
    elpd_pointwise <- loo_obj$pointwise
    
    elpd_pointwise[out_df$ll_ind,1] <- out_df$lpd
    
    return(tibble(country=cntry, 
                  param_type_1=pt1,
                  param_type_2=pt2,
                  looic=-2 * sum(elpd_pointwise[,1]),
                  looic_se=loo_obj$estimates[3,2]))
    
  }
}


# Extract predictive log density for selected points, for a given model index
# corresponding to a line in the lko csv file
get_lpd <- function(model_index, ll_inds,cntry_run_dir){
  file_name <- list.files(cntry_run_dir, 
                          pattern= paste0("*_", model_index,".rds"))
  res <- readRDS(file.path(cntry_run_dir, file_name))
  ll <- extract_log_lik(res$fert_fit)
  
  ll_inds <- as.numeric(str_split(ll_inds, "_")[[1]])
  n_samples <- dim(ll)[1]
  lpd <- apply(ll[,ll_inds, drop=F],2, logSumExp) - log(n_samples)
  return(tibble(ll_ind=ll_inds, lpd=lpd))
}
