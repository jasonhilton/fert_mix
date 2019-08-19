functions{

  // parameteric functions inserted  through R before compilation
  ${param_fun}
}

data {
  int N;
  int n_ages;
  int n_ages_obs;

  
  
  int T_total;
  int T_insample;
  

  int minimum_age;
  int data_gap_low;
  int data_gap_high;
  
  int age_ind[N];
  int time_ind[N];
  
  
  int<lower=0> births[N];
  matrix [n_ages, T_insample] expos;

  int n_keep;
  int keep_ind[n_keep];

  ${dispersion_data}

}

transformed data {
  row_vector[n_ages] ages;

  // We want to evaluate parametric density function at midpoints of each 
  // single-year age group
  // to approximate the mass within each age interval
  // Create a vector of age midpoints to facilitate this.
  for (x in 1:n_ages){
    ages[x] = x - 0.5;
  }

}

parameters {
  vector<lower=0>[T_total] CEB;
  real<lower=0> sigma_CEB;
  
  vector<lower=1>[T_total] sd_age[2];
  real<lower=0> sigma_sd[2];

  // enforce seperation of the modes;
  vector<lower=2>[T_total] mode_age_gap;
  // the minimum of the sum of modes has to be more than the first data age plus
  // the gap to avoid init problems
  vector<lower=36>[T_total] mode_age_sum;
  real<lower=0> sigma_mode[2];
  
  vector[T_total] logit_mixture;
  real<lower=0> sigma_mixture;

  // dispersion parameters inserted before compilation
  ${dispersion_parameters}
}

transformed parameters {
  
  vector[T_total] mode_age_1;
  vector[T_total] mode_age_2;


  vector[N]  log_lik;

  // transformations of dispersion parameters to be inserted before compilation
  ${dispersion_transform}


  mode_age_1 = (mode_age_sum - mode_age_gap) / 2.0;
  mode_age_2 = (mode_age_sum + mode_age_gap) / 2.0;

  {
  matrix[n_ages, T_insample] E_births;
  vector[T_total] mixture;
  matrix[n_ages, T_total] rate;
  ${dispersion_sum_dec}

  for (t in 1:T_total){
    mixture[t] = inv_logit(logit_mixture[t]);
    
    rate[1:n_ages, t] = (CEB[t] * mixture[t] * 
      ${fun_name_1}(ages, mode_age_1[t] - minimum_age, sd_age[1][t]))';

    rate[1:n_ages, t] += (CEB[t] * (1 - mixture[t])* 
      ${fun_name_2}(ages, mode_age_2[t] - minimum_age, sd_age[2][t]))';

  }



  E_births = 1e-06 + rate[1:n_ages, 1:T_insample] .* expos;
  ${dispersion_sums}

   
  for (t in 1:T_insample){
    // avoiding inefficient deep copy through local variable birth_sum
    real birth_sum;
    birth_sum = sum(E_births[1:(data_gap_low + 1),t]);
    E_births[data_gap_low + 1, t] = birth_sum;
    birth_sum = sum(E_births[(n_ages_obs + data_gap_low):n_ages, t]);
    E_births[n_ages_obs + data_gap_low, t] = birth_sum;
  }
  


  for (i in 1:N){
    log_lik[i] = ${lik_dist}(births[i] |
                             E_births[age_ind[i] + data_gap_low, 
                                              time_ind[i]] 
                                      ${dispersion_lik}${dispersion_lik_ind});
  }

  }



}

model {
  // variance priors ------------------------------------------------------------
  
  sigma_CEB ~ normal(0, 4);
  sigma_mode ~ normal(0, 4);
  sigma_sd ~ normal(0, 4);
  sigma_mixture ~ normal(0, 4);

  // dispersion prior inserted before compilation
  ${dispersion_distribution}


  // CEB -----------------------------------------------------------------------
  
  CEB[2:T_total] ~ normal(CEB[1:(T_total - 1)],
                                  sigma_CEB * 0.1);

  CEB[1] ~ normal(2.1, 0.3)T[0,];
  
  // mode ----------------------------------------------------------------------
  mode_age_sum[2:T_total] ~ normal(mode_age_sum[1:(T_total-1)],
                                           sigma_mode[1] * 0.1);

  mode_age_gap[2:T_total] ~ normal(mode_age_gap[1:(T_total-1)],
                                   sigma_mode[2] * 0.1);

  mode_age_sum[1] ~ normal(58,2.5);

  mode_age_gap[1] ~ normal(8,2);  

  // sds -----------------------------------------------------------------------
  for (mm in 1:2){
    sd_age[mm][1]  ~ normal(5, 2)T[1,];
    sd_age[mm][2:T_total] ~ normal(sd_age[mm][1:(T_total - 1)],
                                           sigma_sd[mm] * 0.1);
  }
  
  // logit ---------------------------------------------------------------------

  logit_mixture[2:T_total] ~ normal(logit_mixture[1:(T_total - 1)],
                                            sigma_mixture * 0.1);

  logit_mixture[1] ~ normal(0, 5);

  target += sum(log_lik[keep_ind]);
}
