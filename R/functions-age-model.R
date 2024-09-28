# Functions for AGE FOI model:

read_mumps_data <- function(file){
  # read data and prepare for stan model fitting:
  read.csv(file) %>% mutate(total = positive + negative) %>% 
    mutate(a = 1 + positive, b = 1 + negative,
           low = qbeta(p = c(0.025), a, b), up = qbeta(p = c(0.975), a, b),
           med = qbeta(p = c(0.5), a, b))
           
}

fit_age_model <- function(data_stan){
  age_model <- rstan::stan_model("stan/age_foi_model.stan")
  initfn <- function() { list(log_foi = rep(-8, max(data_stan$chunks))) }
  fit <- rstan::optimizing(age_model, data = data_stan, init = initfn, as_vector = FALSE)
  initf <- function(chain_id = 1) { list( log_foi = fit$par$log_foi ) }
  n_chains <- 4
  init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))
  
  rstan::sampling(age_model, data = data_stan, chains = 4, init = init_ll,
                  iter = 3000, warmup = 900, refresh = 0, seed = 345,
                  control = list(adapt_delta = 0.9999, max_treedepth = 25))
  #return(model_fit)
}

make_fit_summary <- function(model_fit){
  as.data.frame(rstan::summary(model_fit, probs = c(0.025, 0.975))$summary)
  #return(s)
}


# model fit 1:
stan_data_age_modelV1 <- function(mumps_data){
  # model fitting - w/o seroreversion
  list(n_obs = nrow(mumps_data), 
       n_pos = mumps_data$positive,
       n_total = mumps_data$total, 
       ages = mumps_data$Age,
       is_binomial = TRUE, 
       include_seroreversion = 0,
       chunks = mumps_data$Age, 
       age_max = max(mumps_data$Age),
       # these parameters have different meanings dependent on prior choice
       foi_prior_choice = 1, # this prior choice works, w/o errors!
       foi_prior_a = 0, foi_prior_b = 1, 
       serorev_prior_choice = 1, serorev_prior_a = 0, serorev_prior_b = 1
  )
  #returns a list data_stan
}


plot_age_model_fitA <- function(fit_summary, data_stan, mumps_data){
  model_ppc1 <- tibble::rownames_to_column(fit_summary, "parameter") %>% 
    filter(grepl("pos_pred", parameter))
  
  tibble(age = data_stan$ages, actual_prev = (data_stan$n_pos/data_stan$n_total)*100,
         lower = mumps_data$low*100, upper = mumps_data$up*100,
         model_prev = (model_ppc1$mean/data_stan$n_total)*100, 
         ymin = (model_ppc1$`2.5%`/data_stan$n_total)*100, 
         ymax = (model_ppc1$`97.5%`/data_stan$n_total)*100) %>%
    ggplot(aes(x = age)) + 
    geom_line(aes(y = model_prev), color = "#619CFF") + 
    geom_pointrange(aes(y = actual_prev, ymin = lower, ymax = upper), size = 0.2) + 
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#619CFF", alpha = 0.3) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(y = "Seroprevalence (%)", x = "") +
    theme_bw()
  
}

plot_age_model_APIA <- function(mumps_data, fit_summary){
  # annual probability of infection: when FOI is unique for each age:
  ages <- mumps_data$Age
  foi_v1 <- data.frame(age = ages, mean = rep(NA, times = max(ages)), 
                       low = rep(NA, max(ages)), up = rep(NA, max(ages)))
  
  model_foi_v1 <- tibble::rownames_to_column(fit_summary, "parameter") %>% 
    filter(grepl("\\bfoi", parameter))
  
  foi_v1$mean <- 1-exp(-model_foi_v1$mean)
  foi_v1$low <- 1-exp(-model_foi_v1$`2.5%`)
  foi_v1$up <- 1-exp(-model_foi_v1$`97.5%`)
  
  ggplot(foi_v1, aes(x=age)) + 
    geom_line(aes(y = mean), linewidth = 0.7, col = 'black') +
    geom_ribbon(aes(ymin = low, ymax = up), fill = "#1B9E77", alpha = 0.3) +
    labs(x = "", y = "Annual probability\nof infection") + 
    theme_bw()
  
}


# Model 2:
stan_data_age_modelV2 <- function(mumps_data){
  # when FOIs are 'chunks'; and w/o seroreversion:
  ages = mumps_data$Age
  chunks <- rep(c(1,2,3,4), c(2, 8, 10, (max(ages)-(2+8+10)))) # v2 model - chunked FOI
  list(n_obs = nrow(mumps_data),
       n_pos = mumps_data$positive,
       n_total = mumps_data$total,
       ages = ages,
       is_binomial = TRUE,
       include_seroreversion = 0,
       chunks = chunks,
       age_max = max(ages),
       # these parameters have different meanings dependent on prior choice
       foi_prior_choice = 1, # this prior choice works, w/o errors!
       foi_prior_a = 0, foi_prior_b = 1,
       serorev_prior_choice = 1, serorev_prior_a = 0, serorev_prior_b = 1
  )
  
}


plot_age_model_fitB <- function(fit_summary, data_stan, mumps_data){
  model_ppc2 <- tibble::rownames_to_column(fit_summary, "parameter") %>% 
    filter(grepl("pos_pred", parameter))
  
  tibble(age = data_stan$ages, actual_prev = (data_stan$n_pos/data_stan$n_total)*100,
         lower = mumps_data$low*100, upper = mumps_data$up*100,
         model_prev = (model_ppc2$mean/data_stan$n_total)*100, 
         ymin = (model_ppc2$`2.5%`/data_stan$n_total)*100, 
         ymax = (model_ppc2$`97.5%`/data_stan$n_total)*100) %>%
    ggplot(aes(x = age)) + 
    geom_line(aes(y = model_prev), color = "#619CFF") + 
    geom_pointrange(aes(y = actual_prev, ymin = lower, ymax = upper), size = 0.2) + 
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#619CFF", alpha = 0.3) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(y = "Seroprevalence (%)", x = "Age, years") +
    theme_bw()
  
}


plot_age_model_APIB <- function(mumps_data, fit_summary){
  # API: when FOI are chunks for age groups:
  ages <- mumps_data$Age
  foi_v2 <- data.frame(age = ages, mean = rep(NA, times = max(ages)), 
                       low = rep(NA, max(ages)), up = rep(NA, max(ages)))
  
  model_foi_v2 <- tibble::rownames_to_column(fit_summary, "parameter") %>% 
    filter(grepl("\\bfoi", parameter))
  
  foi_v2$mean <- rep(c(model_foi_v2$mean[1],model_foi_v2$mean[2],model_foi_v2$mean[3],model_foi_v2$mean[4]), 
                     c(2, 8, 10, (max(ages)-(2+8+10))))
  foi_v2$low <- rep(c(model_foi_v2$`2.5%`[1],model_foi_v2$`2.5%`[2],model_foi_v2$`2.5%`[3],model_foi_v2$`2.5%`[4]), 
                    c(2, 8, 10, (max(ages)-(2+8+10))))
  foi_v2$up <- rep(c(model_foi_v2$`97.5%`[1],model_foi_v2$`97.5%`[2],model_foi_v2$`97.5%`[3],model_foi_v2$`97.5%`[4]), 
                   c(2, 8, 10, (max(ages)-(2+8+10))))
  
  foi_v2$mean <- 1-exp(-foi_v2$mean)
  foi_v2$low <- 1-exp(-foi_v2$low)
  foi_v2$up <- 1-exp(-foi_v2$up)
  
  ggplot(foi_v2, aes(x=age)) + 
    geom_line(aes(y = mean), linewidth = 0.7, col = 'black') +
    geom_ribbon(aes(ymin = low, ymax = up), fill = "#1B9E77", alpha = 0.3) +
    labs(x = 'Age (years)', y = "Annual probability\nof infection") + 
    theme_bw()
  
}

plots_age_model_fig <- function(A,B,C,D){
  (A + B + C + D) + plot_layout(ncol = 2) + plot_annotation(tag_levels = 'A')
}