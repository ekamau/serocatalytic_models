// constant FOI model: including maternal antibodies

functions{
  array[] real prob_infection_calc(array[] real ages, real lambda, real gamma, int n_obs) {
    array[n_obs] real prob;
    for(i in 1:n_obs){
      real a = ages[i];
      prob[i] = 1-(gamma/gamma-lambda)*(exp(-lambda*a) - exp(-gamma*a));
    }
    
    return prob;
  }
  
}

data{
  int<lower=0> n_obs; // No. rows in data or no. age classes
  array[n_obs] int n_pos; // seropositive
  array[n_obs] int n_total; // tested
  array[n_obs] real ages;
  
}

parameters{
  real<lower=1> gamma; // rate of decay of maternal antibodies
  real<lower=0> lambda; // real<lower=0, upper=gamma>; // 

}

transformed parameters{
  array[n_obs] real<lower=0> prob_infection = prob_infection_calc(ages, lambda, gamma, n_obs);
  array[n_obs] real<lower=0> prob_infectionB;
  
  for(i in 1:n_obs){
    prob_infectionB[i] = 1-exp(-lambda*ages[i]);
  }
    
}

model{
  // priors
  lambda ~ exponential(1);
  gamma ~ cauchy(0, 1);
   
  // likelihood
  n_pos ~ binomial(n_total, prob_infection);
}

generated quantities {
  array[n_obs] int pos_pred = binomial_rng(n_total, prob_infection);
  array[n_obs] int pos_predB = binomial_rng(n_total, prob_infectionB);
  
}
