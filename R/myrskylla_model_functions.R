
# fit the myrskylla et all model to a dataframe containing data from a single country.
fit_single_myrskylla_model <- function(country_data, n_trend=5, n_holdback=10){
  asfr_mat <- country_data %>%  spread(Year, ASFR) %>% 
    select(-Age, -OpenInterval) %>%
    as.matrix()
  
  n_years <- dim(asfr_mat)[2]
  
  first_trend <- n_years - n_holdback - n_trend + 1
  last_trend <- n_years - n_holdback
  trend_inds <- first_trend:last_trend
  ax <- asfr_mat[,last_trend]
  bx <- (asfr_mat[,last_trend] - asfr_mat[,first_trend])/(n_trend - 1)
  A <- asfr_mat[,trend_inds] - ax
  kt <- t(solve(bx %*% bx) %*% bx %*% A)
  epsilon <- diff(kt) - 1
  return(list(ax=ax,
              bx=bx,
              kt=kt,
              epsilon=epsilon,
              delta=1, # by definition
              n_trend=n_trend,
              n_holdback=n_holdback,
              n_years=n_years))
}


# Create predictions from the Myrskylla model, with uncertainty
predict_myrskylla_model <- function(model, data, sigma2, n_predict, n_extrapolate){
  kt <- c(1:n_extrapolate, rep(n_extrapolate, n_predict - n_extrapolate))
  rates <- model$ax  + model$bx %*% t(as.vector(kt))
  
  delta_var <- sigma2 / (model$n_trend - 1)
  extrap_var <- as.matrix(model$bx**2)  %*% (kt**2 * delta_var +  
                                               1:n_predict * sigma2)
  
  qqs <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
  variants <-outer(sqrt(extrap_var), qnorm(qqs))
  for (i in 1:length(qqs)){
    variants[,,i] <- rates + variants[,,i]
  }
  
  jo <- max(data$Year) - n_predict
  Years <- (jo + 1):(jo + n_predict)
  Ages <- unique(data$Age)
  tidy_pred <- cbind(expand.grid(Ages, Years, qqs), as.vector(variants))
  colnames(tidy_pred) <- c("Age", "Year", "Quantile", "ASFR")
  tidy_pred %<>% as_tibble() %>% 
    mutate(Interval=as.factor(abs(Quantile - 0.5) * 2))
  return(tidy_pred)
}
