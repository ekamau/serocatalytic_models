// age-dependent FOI model
// estimates seroreversion parameter

functions{
  real age_foi_calc(int age, int[] chunks, vector foi, real mu) {
    real prob = 0.0;
    for(j in 1:age){
      real lambda = foi[chunks[j]];
      prob = (1 / (lambda + mu)) * exp(-(lambda + mu)) * (lambda * (exp(lambda + mu)  - 1) + prob * (lambda + mu));
    }
    
    return prob;
  }
 
  real[] prob_infection_calc(int[] ages, int[] chunks, vector foi, real mu, int n_obs) {
    real prob_infected[n_obs];
    for(i in 1:n_obs){
      int age = ages[i];
      prob_infected[i] = age_foi_calc(age, chunks, foi, mu);
    }
    
    return prob_infected;
  }
  
}

data{
  int<lower=0> n_obs; // No. rows in data or no. age classes
  int n_pos[n_obs]; // seropositive
  int n_total[n_obs]; // tested
  int<lower=0> age_max;
  int chunks[age_max]; // vector of length len_chunks
  int ages[n_obs];
  
  // model type
  int<lower=0, upper=1> include_seroreversion;
  
  // prior choices
  int<lower=1, upper=6> foi_prior_choice;
  real<lower=0> foi_prior_a; 
  real<lower=0> foi_prior_b;
  
  int<lower=1, upper=3> serorev_prior_choice;
  real<lower=0> serorev_prior_a; 
  real<lower=0> serorev_prior_b;
  
}

transformed data{
  int is_random_walk = foi_prior_choice <= 4 ? 1 : 0;
  int n_chunks = max(chunks); // max value in the vector n_chunks
  
}

parameters{
  row_vector[n_chunks] log_foi; // length of vector = max value in the vector n_chunks
  real<lower=0> sigma[is_random_walk ? 1 : 0]; // only for R/W models
  real<lower=0> nu[foi_prior_choice == 3 ? 1 : 0]; 
  real<lower=0> seroreversion_rate[include_seroreversion ? 1 : 0]; // rate of seroreversion

}

transformed parameters{
  real<lower=0> mu;
  real<lower=0> prob_infection[n_obs];
  vector<lower=0>[n_chunks] foi = to_vector(exp(log_foi));
  
  if(include_seroreversion){
    mu = seroreversion_rate[1];
  } else{
    mu = 0.0;
  }
  
  prob_infection = prob_infection_calc(ages, chunks, foi, mu, n_obs);
    
}

model{
  // likelihood
  n_pos ~ binomial(n_total, prob_infection);
  
  // priors - seroreversion
  if(include_seroreversion) {
    if(serorev_prior_choice == 1) {
      seroreversion_rate ~ cauchy(serorev_prior_a, serorev_prior_b);
      
    } else if (serorev_prior_choice == 2) {
      seroreversion_rate ~ normal(serorev_prior_a, serorev_prior_b);
      
    } else if (serorev_prior_choice == 3) {
      seroreversion_rate ~ uniform(serorev_prior_a, serorev_prior_b);
      
    }
  }
  
  // priors - FOI
  if(foi_prior_choice == 1) { // forward random walk

    sigma ~ cauchy(0, 1);
    log_foi[1] ~ normal(foi_prior_a, foi_prior_b);

    for(i in 2:n_chunks)
      log_foi[i] ~ normal(log_foi[i - 1], sigma);

  } else if(foi_prior_choice == 2) { // backward random walk

    sigma ~ cauchy(0, 1);
    log_foi[n_chunks] ~ normal(foi_prior_a, foi_prior_b);

    for(i in 1:(n_chunks - 1))
      log_foi[n_chunks - i] ~ normal(log_foi[n_chunks - i + 1], sigma);

  } else if(foi_prior_choice == 3){ // forward random walk with Student-t

    sigma ~ cauchy(0, 1);
    nu ~ cauchy(0, 1);
    log_foi[1] ~ normal(foi_prior_a, foi_prior_b);

    for(i in 2:n_chunks)
      log_foi[i] ~ student_t(nu, log_foi[i - 1], sigma);

  } else if (foi_prior_choice == 4) { // backward random walk with Student t

    sigma ~ cauchy(0, 1);
    nu ~ cauchy(0, 1);
    log_foi[n_chunks] ~ normal(foi_prior_a, foi_prior_b);

    for(i in 1:(n_chunks - 1))
      log_foi[n_chunks - i] ~ student_t(nu, log_foi[n_chunks - i + 1], sigma);

  } else if(foi_prior_choice == 5){ // uniform

    foi ~ uniform(foi_prior_a, foi_prior_b);
    target += sum(foi);

  } else if(foi_prior_choice == 6) { // weakly informative

    foi ~ cauchy(foi_prior_a, foi_prior_b);
    target += sum(foi);
   
  }
}

generated quantities {
  int pos_pred[n_obs];
  pos_pred = binomial_rng(n_total, prob_infection);
  
}
