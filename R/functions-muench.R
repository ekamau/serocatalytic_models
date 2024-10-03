
clean_and_name <- function(filename, location) {

  read.csv(filename) %>%
    mutate(
      region=location
      ) %>%
    mutate(
      seropositive=seropositive/100
      )
}

clean_and_name_colombia <- function(filename, location) {
  
  read.csv(filename) %>%
    mutate(
      region=location
    ) %>%
    mutate(seropositive=n_seropositive/n_sample)
}


mle_condition <- function(lambda, Ns, Xs, as) {
  sum((Ns-Xs) * as) - sum(Xs * as * exp(-lambda * as) / (1 - exp(-lambda * as)))
}

estimate_lambda_constant_foi <- function(df) {
  n_sample <- 100 # doesn't matter but can't seem to remove N from condition
  f <- function(lambda) {
    g <- mle_condition(lambda, rep(n_sample, nrow(df)), df$seropositive*n_sample, df$mid_age)
    g^2
  }
  opt <- optim(0.1, f) # only specify initial parameter
  opt$par
}

# assumes single epidemic year
estimate_lambda_colombia <- function(df) {
  - log(1 - sum(df$n_seropositive) / sum(df$n_sample))
}

serocatalytic_constant_foi <- function(age, lambda) {
  1 - exp(-lambda * age)
}

age_sim <- seq(0, 65, 0.1)
create_amazonas_sim <- function(location, lambda_amazonas) {
  
  tibble(mid_age=age_sim) %>%
    mutate(
      modelled=serocatalytic_constant_foi(mid_age, lambda_amazonas)
    ) %>%
    mutate(
      region=location
    )
}

colombian_seroprevalence <- function(age, lambda_colombia) {
  if(age < 5)
    0
  else
    1 - exp(-lambda_colombia)
}

create_colombia_sim <- function(location, lambda_colombia) {
  colombia_sim <- tibble(mid_age=age_sim) %>%
    mutate(
      modelled=map_dbl(mid_age, ~colombian_seroprevalence(., lambda_colombia))
    ) %>%
    mutate(
      region=location
    )
}

combine_all_and_pivot <- function(amazonas_1, amazonas_2, colombia,
                                  amazonas_1_sim, amazonas_2_sim, colombia_sim,
                                  location_1, location_2, location_3) {
  
  df <- amazonas_1 %>%
    bind_rows(
      amazonas_2,
      colombia
    ) %>%
    mutate(
      mid_age=0.5 * (age_lower + age_upper)
    ) %>%
    bind_rows(
      amazonas_1_sim,
      amazonas_2_sim,
      colombia_sim
    ) %>%
    mutate(
      region=as.factor(region)
    ) %>%
    mutate(
      region=fct_relevel(region, location_1, location_2, location_3)
    ) %>%
    mutate(
      shape=1+n_seropositive,
      rate=1+n_sample-n_seropositive
    ) %>%
    mutate(
      lower=qbeta(0.025, shape, rate),
      upper=qbeta(0.975, shape, rate)
    )
  
  df %>%
    select(mid_age, seropositive, modelled, lower, upper, region) %>%
    rename(
      actual=seropositive
    ) %>%
    pivot_longer(c(actual, modelled)) %>%
    mutate(
      name=as.factor(name),
      name=fct_rev(name)
    )
}

plot_muench_seroprevalence <- function(df) {
  ggplot(df, aes(x=mid_age, y=value, group=name)) +
    geom_pointrange(data=df %>% filter(name=="actual"),
                    aes(ymin=lower, ymax=upper)) +
    geom_path(aes(linetype=name)) +
    scale_y_continuous(labels = scales::percent) +
    scale_linetype_manual("Series", breaks=c("actual","modelled"), values=c(1,2)) +
    facet_wrap(~region) +
    theme_classic() +
    theme(
      legend.position = c(0.55, 0.3)
    ) +
    xlab("Age, years") +
    ylab("Seropositive, %")
}

plot_muench_foi <- function(location_1, location_2, location_3,
                            lambda_amazonas_1, lambda_amazonas_2, lambda_colombia) {
  
  year_now <- 1934
  age_sim_short <- seq(0, 65, 1)
  year_start <- year_now - max(age_sim_short)
  years <- rev(seq(year_start, year_now, 1))
  foi_amazonas_1 <- tibble(
    year=years,
    years_ago=age_sim_short,
    foi=lambda_amazonas_1,
    region=location_1
  )
  foi_amazonas_2 <- tibble(
    year=years,
    years_ago=age_sim_short,
    foi=lambda_amazonas_2,
    region=location_2
  )
  foi_colombia <- tibble(
    year=years,
    years_ago=age_sim_short,
    foi=c(0, 0, 0, 0, 0, lambda_colombia, rep(0, length(age_sim_short) - 6)), # assumes study carried out in 1934
    region=location_3
  )
  foi_df <- foi_amazonas_1 %>%
    bind_rows(
      foi_amazonas_2,
      foi_colombia
    ) %>%
    mutate(
      region=fct_relevel(region, location_1, location_2, location_3)
    )
  
  foi_df %>%
    ggplot(aes(x=year, y=foi)) +
    geom_line(linetype=2) +
    scale_x_continuous(breaks = seq(1870, 1930, 10)) +
    scale_y_sqrt() +
    facet_wrap(~region) +
    xlab("Date") +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank()
    ) +
    ylab("FOI, infections per person per year")
}