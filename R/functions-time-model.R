# Functions for TIME FOI model:
read_chikv_data <- function(file){
  df <- read.csv(file) %>% 
    mutate(a = 1 + positive, b = 1 + total - positive) %>% 
    mutate(lower = qbeta(0.025, a, b),
           upper = qbeta(0.975, a, b))
}

stan_data_time_model <- function(df){ # returns a list:
  ages <- df$Age2
  m_exposure <- matrix(nrow = length(ages), ncol = max(ages))
  
  for(i in ages) {
    n_zeros <- max(ages) - i
    n_ones <- i
    m_exposure[which(ages == i), ] <- c(rep(0, n_zeros), rep(1, n_ones))
  }
  
  # No. of FOI parameters to estimate - distinct FOI for each 5 year chunk:
  chunks <- chunks <- rep(seq(1, 11, 1), each=5)
  n_fois_exposed_per_obs <- rowSums(m_exposure)
  foi_index_start_per_obs <- c(1, 1 + cumsum(n_fois_exposed_per_obs))
  foi_index_start_per_obs <- foi_index_start_per_obs[-length(foi_index_start_per_obs)]
  foi_indices <- unlist(map(seq(1, nrow(m_exposure), 1), ~which(m_exposure[., ]==1)))
  
  # data for model fitting w/o seroreversion
  data_stan <- list(
    n_obs=nrow(df),
    n_pos=df$positive,
    n_total=df$total,
    age_max=max(ages),
    observation_exposure_matrix=m_exposure,
    n_fois_exposed_per_obs=n_fois_exposed_per_obs,
    foi_index_start_per_obs=foi_index_start_per_obs,
    include_seroreversion=0,
    n_fois_exposed=sum(n_fois_exposed_per_obs),
    foi_indices=foi_indices,
    chunks=chunks,
    prior_choice=3,
    prior_a=1,
    prior_b=5
  )
  
}

fit_time_model <- function(data_stan){
  model <- rstan::stan_model("stan/time_foi_model.stan")
  init_fn <- function() { list( log_foi = log(rep(0.01, max(data_stan$chunks))) ) }
  fit <- rstan::optimizing(model, data = data_stan, init = init_fn, as_vector = FALSE)
  initf <- function(chain_id = 1) { list( log_foi = fit$par$log_foi ) }
  n_chains <- 4
  init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))
  
  rstan::sampling(model, data=data_stan, init = init_ll, chains = 4, iter = 3000, warmup = 900, 
                      refresh = 0, control = list(adapt_delta = 0.99999, max_treedepth = 25))
  
}

plot_time_model_fit <- function(fit_summary, data_stan, chikv_data){
  model_ppc1 <- tibble::rownames_to_column(fit_summary, "parameter") %>% 
    filter(grepl("pos_pred", parameter))
  
  tibble(age = chikv_data$Age2, actual_prev = (data_stan$n_pos/data_stan$n_total)*100,
         lower = chikv_data$lower*100, upper = chikv_data$upper*100, 
         model_prev = (model_ppc1$mean/data_stan$n_total)*100, 
         ymin = (model_ppc1$`2.5%`/data_stan$n_total)*100, 
         ymax = (model_ppc1$`97.5%`/data_stan$n_total)*100) %>%
    ggplot(aes(x = age)) + geom_line(aes(y = model_prev), color = "#619CFF") + 
    geom_pointrange(aes(y = actual_prev, ymin = lower, ymax = upper), size = 0.2) + 
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#619CFF", alpha = 0.3) +
    labs(y = "Seroprevalence (%)", x = "Age, years") +
    coord_cartesian(ylim = c(0, 100)) +
    theme_bw()
  
}

plot_time_model_API <- function(chikv_data, time_fit_summary){
  # annual probability of infection: when FOI is unique for each year:
  yr_survey <- 2015
  foi <- data.frame(year = rev(yr_survey - seq(0, max(chikv_data$Age2) - 1)), 
                    mean = rep(NA, times = max(chikv_data$Age2)), 
                    low = rep(NA, max(chikv_data$Age2)), up = rep(NA, max(chikv_data$Age2)))
  
  model_foi <- tibble::rownames_to_column(time_fit_summary, "parameter") %>% 
    filter(grepl("fois_by_year", parameter))
  foi$mean <- 1 - exp(-model_foi$mean)
  foi$low <- 1 - exp(-model_foi$`2.5%`)
  foi$up <- 1 - exp(-model_foi$`97.5%`)
  
  ggplot(foi, aes(x=year)) +
    geom_line(aes(y = mean), linewidth = 0.7, col = 'black') +
    geom_ribbon(aes(ymin = low, ymax = up), fill = "#1B9E77", alpha = 0.3) +
    scale_x_continuous(breaks = seq(from = yr_survey, to = min(foi$year), by = -10)) +
    scale_y_continuous(limits = c(0, 0.1)) +
    labs(x = "Year", y = "Annual probability\nof infection") + 
    theme_bw()
  
}

plots_time_model_fig <- function(A,B){
  (A + B) + plot_layout(ncol = 2) + patchwork::plot_annotation(tag_levels = 'A')
}