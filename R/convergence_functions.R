# calculate split rhat statistic (Gelman-Rubin)
calc_rhat <- function(pars, n_par, n_chains=4){
  pars <- pars[,1:n_par]
  n_samples <- dim(pars)[1]
  M <- (2*n_chains)
  N <- n_samples/M
  
  grand_means <- colMeans(pars)
  
  group_means <-matrix(0, n_chains*2, n_par)
  Sm <- matrix(0,M, n_par)
  Bm <- matrix(0,M, n_par)
  for (i in 1:(M)){
    group_means[i,] <- colMeans(pars[1:N + N*(i-1),])
    for (j in 1:N){
      Sm[i,] <- Sm[i,] + ((group_means[i,] -pars[j + N*(i-1),])**2)/(N-1)
    }
    
    Bm[i,] <- (group_means[i,] - grand_means)**2
  }
  
  B <- colSums(Bm)*(N/(M-1))
  
  W <- colSums(Sm)/M
  rhat <- sqrt(((N-1)/N * W + B/N)/W)
  
  return(rhat)
}

# returns the indicies of the chains to throw away for multimodal
# chains
# This a bit ad-hoc, but is only applied in a small number of cases
# for models which perform relatively poorly anyway.
discard_chain_inds <- function(model_obj){
  
  fert_fit <- model_obj$fert_fit
  pars <- as.matrix(fert_fit)
  n_pars <- get_num_parameters(fert_fit)
  pars <- pars[,1:n_pars]
  
  n_chains <- fert_fit@sim$chains
  n_iters_per_chain <- (fert_fit@sim$iter - fert_fit@sim$warmup)/fert_fit@sim$thin
  
  indices <- lapply(0:3, function(i) 1:n_iters_per_chain +  i * n_iters_per_chain)
  # calc_rhats throwing away one chain at a time
  rhats <- lapply(indices, function(ind){
    calc_rhat(pars[-ind,],n_pars, n_chains=(n_chains-1)*2)
  })
  # which chain, when discarded, results in convergence?
  which_conv <- lapply(rhats, function(rh) sum(rh>1.05)==0)
  
  
  ## what if none of them do?
  if(sum(unlist(which_conv))==0){
    print("More than one chain at different mode...")
    
    # split into groups of two chains and repeat
    ind1 <- c(unlist(indices[[1]]),unlist(indices[[2]]))
    ind2 <- c(unlist(indices[[1]]),unlist(indices[[3]]))
    ind3 <- c(unlist(indices[[1]]),unlist(indices[[4]]))
    indices_2 <- list(ind1, ind2, ind3)
    rhats <- lapply(indices_2, function(ind){
      calc_rhat(pars[-ind,],n_pars, n_chains=(n_chains-1)*2)
    })
    # are the other two  chains in each case converged (to each other)
    rhats_neg <- lapply(indices_2, function(ind){
      calc_rhat(pars[ind,],n_pars, n_chains=(n_chains-1)*2)
    })
    # which of the inital groups are converged
    which_conv <- lapply(rhats, function(rh) sum(rh>1.05)==0)
    which_conv_neg <- lapply(rhats_neg, function(rh) sum(rh>1.05)==0)
    # if nothing is converged, give up
    if((sum(unlist(which_conv))==0) & (sum(unlist(which_conv_neg))==0)){
      return(NA)
    } else if(sum(unlist(which_conv))==0) {
      # if removing only chains from the negative list converge
      # return those
      inds <- unlist(indices_2[unlist(which_conv_neg)])
      all_inds <- (1:(n_iters_per_chain*n_chains))
      return(all_inds[!(all_inds %in% inds)])
    }
    
    if(which_conv_neg[which(unlist(which_conv))[1]]==T){
      # both are converged
      ll <- extract_log_lik(fert_fit, merge_chains = F)
      inds <- unlist(indices_2[unlist(which_conv)])
      loo_select <- compare_by_loos(chains_from_inds(inds, 
                                         n_iters_per_chain * 4), ll)
      
      if (loo_select==T){
        # chains defined by ind have a better looic than the inverse
        # so we want to throw away the inverse
        all_inds <- (1:(n_iters_per_chain*n_chains))
        return(all_inds[!(all_inds %in% inds)])
      } else {
        # these chains perform worse, so throw them away.
        return(inds)
      }
    } else {
      # only removing the first pair results in convergence
      return(unlist(indices_2[unlist(which_conv)]))
    }
  }
  return(unlist(indices[unlist(which_conv)]))
}

# compare the relative performance of different chains
compare_by_loos <- function(chains,ll){
  ll1 <- ll[,chains,]
  r_eff1 <- relative_eff(ll1)
  loo1 <- loo(ll1, r_eff = r_eff1)
  ll2 <- ll[,-chains,]
  r_eff2 <- relative_eff(ll2)
  loo2 <- loo(ll2, r_eff = r_eff2)
  comp <- compare(loo1, loo2)
  # return true if chains are better than not chains
  return(comp[1]<0)
}

# return chain index given iteration index
chains_from_inds <- function(inds,n_iter, n_chains=4){
  iter_p_chain <- floor(n_iter/n_chains)
  return(unique(floor((inds-1) / iter_p_chain))+1)
}

