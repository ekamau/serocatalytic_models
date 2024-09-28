// constant FOI model: including maternal antibodies

functions{
  real[] prob_infection_calc(real[] ages, real lambda, real gamma, int n_obs) {
    real prob[n_obs];
    for(i in 1:n_obs){
      real a = ages[i];
      prob[i] = 1-(gamma/gamma-lambda)*(exp(-lambda*a) - exp(-gamma*a));
    }
    
    return prob;
  }
  
}

data{
  int<lower=0> n_obs; // No. rows in data or no. age classes
  int n_pos[n_obs]; // seropositive
  int n_total[n_obs]; // tested
  real ages[n_obs];
  
}

parameters{
  real<lower=1> gamma; // rate of decay of maternal antibodies
  real<lower=0> lambda; // real<lower=0, upper=gamma>; // 

}

transformed parameters{
  real<lower=0> prob_infection[n_obs] = prob_infection_calc(ages, lambda, gamma, n_obs);
  real<lower=0> prob_infectionB[n_obs];
  
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
  int pos_pred[n_obs] = binomial_rng(n_total, prob_infection);
  int pos_predB[n_obs] = binomial_rng(n_total, prob_infectionB);
  
}
