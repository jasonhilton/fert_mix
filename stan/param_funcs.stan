
  row_vector gamma_density(row_vector xx, real location, real scale){
    // location is the mode
    // scale is the sd.
    real rate;
    real shape;
    int n_ages = num_elements(xx);
    row_vector[n_ages] out;


    rate = (location + sqrt(location^2 + 4 * scale^2)) / (2 * scale^2);
    shape = location * rate + 1;
    for (x in 1:n_ages){
      out[x] = exp(gamma_lpdf(xx[x] | shape, rate)); // allows non integer xx;
    }
    return(out);
  }

  row_vector weibull_density(row_vector xx, real location, real scale){
   // location is the median
   // scale is the distance from the median to the upper quartile. 
  
    real upper_quartile;
    real k;
    real lam;
    int n_ages = num_elements(xx);
    row_vector[n_ages] out;
    upper_quartile = location + scale;
    
     
    k = log(2)/(log(upper_quartile) - log(location));
  
    lam = location / ((log(2))^(1.0/k));
    for (x in 1:n_ages){
      out[x] = exp(weibull_lpdf(xx[x] | k, lam));
    }
   return (out);

}

row_vector hadwiger_density(row_vector xx, real location, real scale){
   // location is the mean
   // scale is the sd
  int n_ages = num_elements(xx);
  row_vector[n_ages] out;
  real part1;
  real part2;
  real part3;
  real log_H;
  
  log_H = log(location) - 0.5*(log(2)) - log(scale);
  part1 = log_H - log(location) - 0.5 * log(pi());

  for (x in 1:n_ages){
    part2 = 1.5  * log(location) - 1.5 * log(xx[x]);
    part3 = - (exp(log_H)^2)  * (location/xx[x] + xx[x]/location - 2);
    out[x] = (exp(part1 + part2 + part3));
  }
  return (out);

}


