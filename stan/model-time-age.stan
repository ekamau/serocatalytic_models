// without seroreversion assumption!
// estimated age rates; fixed time rates; 
// works well with optimizing

functions {
  matrix prepare_mat_t(){
    matrix[12,12] mat_t = rep_matrix(0, 12, 12);
    
    for(i in 1:12){
      mat_t[i,i] = -1;
      if(i > 1){
        mat_t[i,(i-1)] = 1;
      }
    }
    
    mat_t[1,1] = 0; mat_t[2,1] = 0; mat_t[12,12] = 0;
    return(mat_t);
  }
  
  real prob_seropos_calculate(vector time_rate, vector age_rate, int age){
    
    vector[age] v_vector = tail(time_rate, age);
    matrix[12,12] mat_A = rep_matrix(0, 12, 12);
    vector[12] y_vec_init = [1,0,0,0,0,0,0,0,0,0,0,0]';
      
    for(i in 1:age){
      matrix[12,12] mat_t = prepare_mat_t();
      real foi = (age_rate[i] * v_vector[i]); // foi
      mat_t[1,1] = -(foi); 
      mat_t[2,1] = foi; 
      // add matrix to another matrix:
      mat_A += mat_t; // sum of matrices
    }
    
    vector[12] mat_yvec_t = matrix_exp(mat_A) * y_vec_init; # matrix and vector product
    real prob_seropos = sum(mat_yvec_t[2:11]) / (1-mat_yvec_t[12]);
      
    return prob_seropos;
  }
  
  vector probs_generate(vector time_rate, vector age_rate, int N, int[] ages, int[] time_chunks, int age_max){
    
    vector[N] prob;
    vector[age_max] time_rate_longer;
    
    for(i in 1:age_max){
      time_rate_longer[i] = time_rate[time_chunks[i]];
    }
    
    for(n in 1:N){
      int age = ages[n];
      prob[n] = prob_seropos_calculate(time_rate_longer, age_rate, age);
    }
    return prob;
  }
  
}
      
data {
  int<lower=0> N;
  int<lower=0> age_max;
  int n_chunks;
  int time_chunks[age_max];
  int ages[N];
  int n_pos[N];
  int total[N];
 
}

parameters {
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> c;
  real<lower=0> sigma;
  vector<lower=0, upper=1>[n_chunks] time_rate;
  
}

transformed parameters {
  vector<lower=0>[age_max] age_rate;
  
  for(i in 1:age_max){
    age_rate[i] = c*exp(gamma_lpdf(i | a, b));
  }
    
  vector<lower=0>[N] probability_age = probs_generate(time_rate, age_rate, N, ages, time_chunks, age_max);
  
}

model {
  // priors
  a ~ cauchy(0,1);
  b ~ cauchy(0,1);
  c ~ cauchy(0,1);
  sigma ~ cauchy(0,1);
  
  // likelihood
  n_pos ~ binomial(total, probability_age);
  
}

generated quantities {
  int pos_pred[N] = binomial_rng(total, probability_age);
  
}
