# Functions for maternal antibodies model:
read_ev68_data <- function(file){
  # read data and prepare for stan model fitting:
  df <- read.csv("EV68.csv") %>% 
    dplyr::filter(Year == 2006) %>% 
    mutate(seroStatus = case_when(final_Titre >= 16 ~ 'Positive',
                                  final_Titre < 16 ~ 'Negative'),
           age_rounded = case_when(Age < 0.4 ~ 0.2,
                                   (Age >= 0.4) & (Age < 1) ~ 0.7,
                                   TRUE ~ floor(Age)))
  
  df1 <- df %>% count(age_rounded, seroStatus) %>% 
    reshape2::dcast(age_rounded ~ seroStatus, value.var = "n") %>% 
    replace(is.na(.), 0) %>% 
    mutate(total = Negative + Positive,
           prop = Positive/total) %>% 
    mutate(a = 1 + Positive, b = 1 + total - Positive) %>% 
    mutate(lower = qbeta(0.025, a, b),
           upper = qbeta(0.975, a, b))
  
}


stan_data_ev68_model <- function(df){
  list(n_obs = nrow(df),
       n_pos = df$Positive,
       n_total = df$total,
       ages = df$age_rounded
  )
  #returns a list data_stan
  
}


fit_ev68_model <- function(data_stan){
  model <- stan_model("stan/maternal_abs_constant_FOI.stan")
  initfn <- function() { list(lambda = 0.05, gamma = 2) }
  fit <- optimizing(model, data = data_stan, init = initfn, as_vector = FALSE)
  
  initf <- function(chain_id = 1) { list( lambda = fit$par$lambda, gamma = fit$par$gamma ) }
  n_chains <- 4
  init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))
  
  rstan::sampling(model, data = data_stan, chains = 4, init = init_ll, iter = 3000, warmup = 900, 
                  refresh = 0, seed = 345, control = list(adapt_delta = 0.9999, max_treedepth = 25))
  
}


plot_ev68_model_fit <- function(ev68_model_fit, ev68_data, ev68_data_stan){
  s <- as.data.frame(summary(ev68_model_fit, probs = c(0.025, 0.975))$summary)
  model_ppc <- tibble::rownames_to_column(s, "parameter") %>% 
    dplyr::filter(grepl("\\bprob_infection\\b", parameter))
  
  fig1 <- tibble(age=ev68_data_stan$ages, lower=ev68_data$lower*100, upper=ev68_data$upper*100,
                 actual=ev68_data$prop*100, model_predicted=model_ppc$mean*100,
                 ymin=model_ppc$`2.5%`*100, ymax=model_ppc$`97.5%`*100) %>% 
    ggplot(aes(x=age)) +
    geom_pointrange(aes(y=actual, ymin=lower, ymax = upper)) +
    geom_line(aes(y=model_predicted), color = "#619CFF") +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), fill = "#619CFF", alpha = 0.3) +
    coord_cartesian(xlim = c(0, 30), ylim = c(0, 100)) +
    scale_x_sqrt() +
    labs(x = "Age, years", y = "Seroprevalence (%)",
         subtitle = "Considering maternal antibodies") +
    theme_bw()
  
  ## compare fitting with model with constant FOI model w/o considering maternal antibodies
  model_ppc2 <- tibble::rownames_to_column(s, "parameter") %>%
    filter(grepl("\\bprob_infectionB\\b", parameter))

  fig2 <- tibble(age=ev68_data_stan$ages, actual=ev68_data$prop*100, lower=ev68_data$lower*100,
                 upper=ev68_data$upper*100, model_predicted=model_ppc2$mean*100,
                 ymin=model_ppc2$`2.5%`*100, ymax=model_ppc2$`97.5%`*100) %>%
    ggplot(aes(x=age)) + 
    geom_pointrange(aes(y=actual, ymin=lower, ymax = upper)) +
    geom_line(aes(y=model_predicted), color = "#159090") +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), fill = "#159090", alpha = 0.3) +
    coord_cartesian(xlim = c(0, 30), ylim = c(0, 100)) +
    scale_x_sqrt() +
    labs(x = "Age, years", y = "",
         subtitle = "Neglecting maternal antibodies") +
    theme_bw()

  fig1 + fig2 + patchwork::plot_annotation(tag_levels = "A")
  
}