# Functions for combined time and age model - HIV model:

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE, threads_per_chain = 2)

read_hiv_data <- function(file){
  df <- read.csv(file, header = TRUE) %>% 
    select(Age, prev2003, n2003) %>% 
    mutate(prev = prev2003 / n2003)
}


calculate_binomial_ci <- function(df){
  seroprev <- data.frame()
  for(i in 1:nrow(df)){
    ci <- qbeta(p = c(0.025, 0.5, 0.975), 1+(df$prev2003[i]), 1+(df$n2003[i]-df$prev2003[i]))
    seroprev <- rbind(seroprev, c(df$Age[i], df$prev[i], ci))
  }
  
  seroprev <- seroprev %>% `colnames<-`(c("age", "prev", "lower", "median", "upper"))
}


stan_data_time_age_model <- function(df){
  data_stan <- list(
    N=nrow(df),
    age_max=max(df$Age),
    ages=df$Age,
    n_pos=df$prev2003,
    total=df$n2003,
    n_chunks=2,
    time_chunks=rep(c(1,2), c(27,22)))
  
}

### NOT USED / RUN:
# fit_time_age_model <- function(data_stan){
#   u <- 2 * dlnorm(1:max(data_stan$ages), meanlog = 3.5, sdlog = 0.5)
#   model <- stan_model("stan/time_age_foi_model.stan")
#   init_fn <- function() { list(age_rate = u, a = 20, b = 1, c = 1) }
#   fit_optim <- optimizing(model, data = data_stan, init = init_fn, as_vector = FALSE)
#   
#   initf <- function(chain_id = 1) {
#     list( a = fit_optim$par$a, b = fit_optim$par$b, c = fit_optim$par$c,
#           time_rate = fit_optim$par$time_rate, sigma = fit_optim$par$sigma )
#   }
# 
#   n_chains <- 4
#   init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))
# 
#   fit_sample <- sampling(model, data = data_stan, init = init_ll, chains = 4, iter = 3000, warmup = 900,
#                          refresh = 1000, control = list(adapt_delta = 0.99999, max_treedepth = 25))
#   
# }


# panel A of figure in manuscript:
plot_model_fit_hiv <- function(model_fit_hiv, df, data_stan){
  s <- as.data.frame(summary(model_fit_hiv, probs = c(0.025, 0.975))$summary)
  pred_seropos <- tibble::rownames_to_column(s, "parameter") %>% filter(grepl("pos_pred", parameter))
  
  tibble(age = df$age, actual_seroprev = df$prev, model_seroprev = pred_seropos$mean/data_stan$total,
         actual_seroprev_low = df$lower, actual_seroprev_high = df$upper, 
         model_seroprev_low = pred_seropos$`2.5%`/data_stan$total, 
         model_seroprev_high = pred_seropos$`97.5%`/data_stan$total) %>%
    ggplot(aes(x = age)) +
    geom_line(aes(y = model_seroprev), color = "#619CFF", linewidth = 1) + 
    geom_pointrange(aes(y = actual_seroprev, ymin = actual_seroprev_low, ymax = actual_seroprev_high)) +
    geom_ribbon(aes(y = model_seroprev, ymin = model_seroprev_low, ymax = model_seroprev_high), 
                alpha = 0.2, fill = "#619CFF") +
    scale_x_continuous(limits = c(10,50), breaks = seq(0, 50, by=5)) +
    scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
    labs(y = "Seroprevalence", x = "Age, years") + 
    theme_bw()
  
}


# function to construct A matrix for one piece:
prep_matrix <- function(){
  mat <- matrix(0, nrow = 12, ncol = 12) # array(0L, dim=c(12,12))
  for(i in 1:nrow(mat)){ # i represents compartment
    mat[i,i] <- -1
    if(i > 1){
      mat[i,(i-1)] <- 1
    }
  }
  mat[1,1] <- 0; mat[2,1] <- 0; mat[12,12] <- 0
  return(mat)
}


calculate_prob_infection <- function(model_fit_hiv, df){
  # Model estimated time and age rates:
  s <- as.data.frame(summary(model_fit_hiv, probs = c(0.025, 0.975))$summary)
  u_df <- tibble::rownames_to_column(s, "parameter") %>% filter(grepl("age_rate", parameter))
  v_df <- tibble::rownames_to_column(s, "parameter") %>% filter(grepl("time_rate", parameter))
  u <- u_df$mean; v <- v_df$mean
  
  year_survey <- 2003
  ages <- df$Age # original HIV data
  yrs <- seq((year_survey-max(ages)+1), year_survey, 1)
  foi_df_time <- data.frame(year = yrs, foi = c(rep(v[1], (which(yrs == 1980))), 
                                                rep(v[2], (length(yrs) - (which(yrs == 1980))))))
  
  #- ODE solution with matrix exponentials => compute probability of being seropositive
  y_vec_init = c(1,rep(0,11))
  v <- foi_df_time$foi
  prob_seropos_df <- data.frame()
  prob_infection_df <- data.frame()
  
  for(age in ages){
    years = rev(head(year_survey:(year_survey - age), -1)) # years of exposure by the cohort - should be length of age
    cohort = (year_survey - age) # birth year
    v_vector = tail(v, age)
    mat_A = matrix(0, nrow = 12, ncol = 12)
    
    for(t in years){
      mat_t <- prep_matrix() # matrix A at time/year t
      a = t - cohort # age in year t
      foi = (u[a] * v_vector[a])
      mat_t[1,1] = -(u[a] * v_vector[a]) # -foi
      mat_t[2,1] = (u[a] * v_vector[a]) # foi
      # add matrix to another matrix:
      mat_A <- mat_A + mat_t # sum of matrices
      
      prob_infection_df = rbind(prob_infection_df, c(age, cohort, t, foi))
    }
    
    mat_yvec_t = Matrix::expm(mat_A) %*% y_vec_init # matrix and vector product
    prob_seropos = sum(mat_yvec_t[2:11]) / (1-mat_yvec_t[12])
    prob_seropos_df = rbind(prob_seropos_df, c(age, cohort, prob_seropos))
    
  }
  
  prob_infection_df <- prob_infection_df %>% `colnames<-`(c("age","birth-yr","year","foi"))
  prob_seropos_df <- prob_seropos_df %>% `colnames<-`(c("age","birth-yr","prob_seropos"))
  
  return(prob_infection_df)
  
}


# panel C of figure in manuscript:
plot_prob_infection <- function(df){ # prop_infection_df
  df %>% ggplot(aes(x = year, y = (1-exp(-foi)), group = as.factor(age))) +
    geom_line(aes(color = age)) +
    scale_color_continuous(type = "viridis", option = 'A', direction = -1, 
                           limits = c(min(unique(df$age)), max(unique(df$age))),
                           oob = scales::squish, 
                           breaks = seq(min(unique(df$age)), max(unique(df$age)), 5),
                           labels = ~ ifelse(.x < max(unique(df$age)), .x, '49'),
                           guide = guide_colorsteps(show.limits = TRUE)) +
    labs(x = "Year", y = "Prob(infected)", color = "Age, years\n") + 
    theme_bw() +
    theme(legend.text = element_text(size = 8),
          legend.title = element_text(size = 10),
          legend.position = c(0.2,0.6))
  
}


# panel B of figure in manuscript:
plot_age_rate <- function(model_fit_hiv, df){
  s <- as.data.frame(summary(model_fit_hiv, probs = c(0.025, 0.975))$summary)
  u_df <- tibble::rownames_to_column(s, "parameter") %>% 
    filter(grepl("age_rate", parameter))
  foi_df_age <- data.frame(age = 1:max(df$Age), foi = u_df$mean)
  
  ggplot(foi_df_age, aes(x = age, y = foi)) + 
    geom_line() +
    labs(x = "Age, years", y = "Age rate (u)") +
    theme_bw()
  
}


# panel D of figure in manuscript:
plot_time_rate <- function(model_fit_hiv, df){
  s <- as.data.frame(summary(model_fit_hiv, probs = c(0.025, 0.975))$summary)
  v_df <- tibble::rownames_to_column(s, "parameter") %>% filter(grepl("time_rate", parameter))
  year_survey <- 2003
  yrs <- seq((year_survey-max(df$Age)+1), year_survey, 1)
  v <- v_df$mean
  
  data.frame(year = yrs, foi = c(rep(v[1], (which(yrs == 1980))),
                                 rep(v[2], (length(yrs) - (which(yrs == 1980)))))) %>% 
    ggplot(aes(x = year, y = foi)) + 
    geom_line() +
    labs(x = "Year", y = "Time rate (v)") +
    theme_bw()
  
}


# Plot final figure:
plot_model_fit_hiv_fig <- function(figA,figB,figD,figC){
  (figA + figB + plot_layout(widths = c(2,1))) / (figD + figC + plot_layout(widths = c(2,1))) + 
    plot_annotation(tag_levels = "A")
  
}

