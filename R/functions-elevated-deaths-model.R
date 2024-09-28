# Functions for elevated deaths (Ebola) model:

read_ebola_data <- function(file){
  # read data and prepare for stan model fitting:
  read.csv(file) %>% mutate(age = floor((age_min + age_max)/2))
  
}

calculate_binom_int <- function(df){
  df$prev <- df$seropos/df$n
  seroprev <- data.frame()
  
  for(i in 1:nrow(df)){
    ci <- qbeta(p = c(0.025, 0.5, 0.975), 1+df$seropos[i], 1+(df$n[i]-df$seropos[i]))
    seroprev <- rbind(seroprev, c(df$age[i], df$prev[i], ci))
  }
  seroprev <- seroprev %>% `colnames<-`(c("age", "prev", "lower", "median", "upper"))
  
}

stan_data_ebola_model <- function(df){
  list(n_obs = nrow(df),
       n_pos = df$seropos,
       n_total = df$n,
       ages = df$age
  )
  #returns a list data_stan
  
}

fit_ebola_modelA <- function(data_stan){
  model <- rstan::stan_model("stan/elevated_death_modelA.stan")
  initfn <- function() { list(lambda = 0.1, mu = 0.04, rho = 0.89) }
  fit <- optimizing(model, data = data_stan, init = initfn, as_vector = FALSE)
  
  initf <- function(chain_id = 1) { 
    list( lambda = fit$par$lambda, mu = fit$par$mu, rho = fit$par$rho ) }
  n_chains <- 4
  init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))
  
  rstan::sampling(model, data = data_stan, chains = 4, init = init_ll, iter = 3000, warmup = 900, 
                  refresh = 0, seed = 345, control = list(adapt_delta = 0.9999, max_treedepth = 25))
  #return(model_fit)
  
}

fit_ebola_modelB <- function(data_stan){
  model <- rstan::stan_model("stan/elevated_death_modelB.stan")
  initfn <- function() { list(lambda = 0.1, mu = 0.04, rho = 0.0) }
  fit <- optimizing(model, data = data_stan, init = initfn, as_vector = FALSE)
  
  initf <- function(chain_id = 1) { 
    list( lambda = fit$par$lambda, mu = fit$par$mu, rho = fit$par$rho ) }
  n_chains <- 4
  init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))
  
  rstan::sampling(model, data = data_stan, chains = 4, init = init_ll, iter = 3000, warmup = 900, 
                  refresh = 0, seed = 345, control = list(adapt_delta = 0.9999, max_treedepth = 25))
  
}


plot_ebola_model_fit <- function(fitA, fitB, data_stan, seroprev){
  sA <- as.data.frame(summary(fitA, probs = c(0.025, 0.975))$summary)
  ppcA <- tibble::rownames_to_column(sA, "parameter") %>% dplyr::filter(grepl("\\bpos_pred\\b", parameter))
  
  sB <- as.data.frame(summary(fitB, probs = c(0.025, 0.975))$summary)
  ppcB <- tibble::rownames_to_column(sB, "parameter") %>% dplyr::filter(grepl("\\bpos_pred\\b", parameter))

  modA <- tibble(age=data_stan$ages, actual=seroprev$prev, actual_low=seroprev$lower,
                 actual_up=seroprev$upper, model_predicted=ppcA$mean/data_stan$n_total,
                 ymin=ppcA$`2.5%`/data_stan$n_total, ymax=ppcA$`97.5%`/data_stan$n_total) %>%
    mutate(model = "m1")
  
  modB <- tibble(age=data_stan$ages, actual=seroprev$prev, actual_low=seroprev$lower,
                 actual_up=seroprev$upper, model_predicted=ppcB$mean/data_stan$n_total,
                 ymin=ppcB$`2.5%`/data_stan$n_total, ymax=ppcB$`97.5%`/data_stan$n_total) %>%
    mutate(model = "m2")

  rbind(modA, modB) %>% ggplot(aes(x=age)) +
    geom_pointrange(aes(y = actual, ymin = actual_low, ymax = actual_up)) +
    geom_line(aes(y = model_predicted, color = model), position=position_jitter(height=0.005)) +
    scale_color_manual(labels = c("IFR = 0.89", "IFR = 0.0"),
                       values = c("m1" = "#619CFF", "m2" = "#159090")) +
    scale_y_continuous(labels = function(x) x*100, limits = c(0,1)) +
    labs(x = "Age, years", y = "Seroprevalence (%)", color = "") +
    theme_bw() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.7,0.8))
  
}

plot_ebola_model_API <- function(fitA, fitB){
  sA <- as.data.frame(summary(fitA, probs = c(0.025, 0.975))$summary)
  sB <- as.data.frame(summary(fitB, probs = c(0.025, 0.975))$summary)
  
  # annual probability of infection: 
  lambdaA <- (tibble::rownames_to_column(sA, "parameter") %>% filter(grepl("\\blambda\\b", parameter)))
  lambdaB <- (tibble::rownames_to_column(sB, "parameter") %>% filter(grepl("\\blambda\\b", parameter)))
  
  bind_rows((lambdaA %>% mutate(model = "A")), (lambdaB %>% mutate(model = "B"))) %>% 
    mutate(mean = 1-exp(-mean), `2.5%` = 1-exp(-`2.5%`), `97.5%` = 1-exp(-`97.5%`)) %>% 
    ggplot(aes(x = model, y = mean, color = model)) +
    geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
    scale_color_manual(breaks = c("A", "B"), values = c("#619CFF","#159090")) +
    scale_x_discrete(name ="Model", labels = c("A"="IFR = 0.89", "B"="IFR = 0.0")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.2)) +
    ylab("Annual probability\nof infection") +
    theme_bw() +
    theme(legend.position = "none")
  
}

plots_ebola_model_fig <- function(A,B){
  (A + B) + plot_annotation(tag_levels = 'A')
}