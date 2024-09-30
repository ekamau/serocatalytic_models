# Created by use_targets().
# Load packages required to define the pipeline:
library(targets)
library(crew)

# Set target options:
tar_option_set(
  packages = c("tibble","tidyverse","rstan","patchwork","Hmisc","Matrix","expm"), 
  controller = crew::crew_controller_local(workers = 4, seconds_idle = 60)

  # Pipelines that take a long time to run may benefit from optional distributed computing. 
  # To use this capability in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  # 
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to: workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
)

# Run the R scripts in the R/ folder with custom functions:
tar_source("R/functions.R")

# Needed for naming Muench-type analyses
name_1 <- "i. Taracua, Sao Gabriel and Yaurete\nBrazil"
name_2 <- "ii. Esperanga and Regiao do Sul de Tabatinga\nBrazil"
name_3 <- "iii. Socorro\nColombia"

# target list:
list(
  ### Age FOI model:
  #1) estimate FOI for each age unit:
  tar_target(mumps_data, read_mumps_data(file = "data/mumps.csv")),
  tar_target(mumps_data_stanV1, stan_data_age_modelV1(mumps_data)),
  tar_target(age_model_fitV1, fit_age_model(mumps_data_stanV1)),
  tar_target(age_fit_summaryV1, make_fit_summary(age_model_fitV1)),
  tar_target(plot_age_fitA, plot_age_model_fitA(age_fit_summaryV1, mumps_data_stanV1,
                                               mumps_data)),
  tar_target(plot_age_apiA, plot_age_model_APIA(mumps_data, age_fit_summaryV1)),

  #2) estimate FOI for age groups:
  tar_target(mumps_data_stanV2, stan_data_age_modelV2(mumps_data)),
  tar_target(age_model_fitV2, fit_age_model(mumps_data_stanV2)),
  tar_target(age_fit_summaryV2, make_fit_summary(age_model_fitV2)),
  tar_target(plot_age_fitB, plot_age_model_fitB(age_fit_summaryV2, mumps_data_stanV2,
                                               mumps_data)),
  tar_target(plot_age_apiB, plot_age_model_APIB(mumps_data, age_fit_summaryV2)),
  tar_target(plot_age_FOI_model,
             plots_age_model_fig(plot_age_fitA, plot_age_apiA, plot_age_fitB, plot_age_apiB)),
  
  ### Time FOI model:
  tar_target(chikv_data, read_chikv_data(file = "data/chikv.csv")),
  tar_target(chikv_data_stan, stan_data_time_model(chikv_data)),
  tar_target(time_model_fit, fit_time_model(chikv_data_stan)),
  tar_target(time_fit_summary, make_fit_summary(time_model_fit)),
  tar_target(plot_time_fit, plot_time_model_fit(time_fit_summary, chikv_data_stan, chikv_data)),
  tar_target(plot_time_api, plot_time_model_API(chikv_data, time_fit_summary)),
  tar_target(plot_time_FOI_model, plots_time_model_fig(plot_time_fit, plot_time_api)),

  ### Elevated deaths model:
  tar_target(ebov_data, read_ebola_data(file = "data/ebola.csv")),
  tar_target(ebov_data_stan, stan_data_ebola_model(ebov_data)),
  tar_target(ebola_model_fitA, fit_ebola_modelA(ebov_data_stan)),
  tar_target(ebola_model_fitB, fit_ebola_modelB(ebov_data_stan)),
  tar_target(seroprev_data, calculate_binom_int(ebov_data)),
  tar_target(plot_ebola_fit, plot_ebola_model_fit(ebola_model_fitA, ebola_model_fitB,
                                                  ebov_data_stan, seroprev_data)),
  tar_target(plot_ebola_api, plot_ebola_model_API(ebola_model_fitA, ebola_model_fitB)),
  tar_target(plot_ebola_model_result, plots_ebola_model_fig(plot_ebola_fit, plot_ebola_api)),
  
  ### maternal antibodies FOI model:
  tar_target(ev68_data, read_ev68_data(file = "data/EV68.csv")),
  tar_target(ev68_data_stan, stan_data_ev68_model(ev68_data)),
  tar_target(ev68_model_fit, fit_ev68_model(ev68_data_stan)),
  tar_target(plot_ev68_fit, plot_ev68_model_fit(ev68_model_fit, ev68_data, ev68_data_stan)),
  
  ### Time and Age model (HIV model):
  tar_target(hiv_data, read_hiv_data(file = "data/mossong_HIV.csv")),
  tar_target(hiv_data2, calculate_binomial_ci(hiv_data)),
  tar_target(hiv_data_stan, stan_data_time_age_model(hiv_data)),
  #** files uploaded to Github don't have targets below: **
  tar_target(hiv_model_fit, fit_time_age_model(hiv_data_stan)),
  tar_target(plot_hiv_fit, plot_model_fit_hiv(hiv_model_fit, hiv_data, hiv_data_stan)),
  tar_target(hiv_prob_infxn, calculate_prob_infection(hiv_model_fit, hiv_data)),
  tar_target(plot_hiv_prob_infxn, plot_prob_infection(hiv_prob_infxn)),
  tar_target(plot_age_foi, plot_age_rate(hiv_model_fit, hiv_data)),
  tar_target(plot_time_foi, plot_time_rate(hiv_model_fit, hiv_data)),
  tar_target(plot_hiv_model_fig, plot_model_fit_hiv_fig(plot_hiv_fit, plot_age_foi,
                                                        plot_hiv_prob_infxn, plot_time_foi)),
  
  # Muench data analysis
  
  ## process raw data
  tar_target(df_amazonas_1, clean_and_name("data/amazonas_1.csv", name_1)),
  tar_target(df_amazonas_2, clean_and_name("data/amazonas_2.csv", name_2)),
  tar_target(df_colombia, clean_and_name_colombia("data/colombia.csv", name_3)),
  
  ## estimate lambdas
  tar_target(lambda_amazonas_1, estimate_lambda_constant_foi(df_amazonas_1)),
  tar_target(lambda_amazonas_2, estimate_lambda_constant_foi(df_amazonas_2)),
  tar_target(lambda_colombia, estimate_lambda_colombia(df_colombia)),
  
  ## process to plot
  tar_target(df_amazonas_sim_1, create_amazonas_sim(name_1, lambda_amazonas_1)),
  tar_target(df_amazonas_sim_2, create_amazonas_sim(name_2, lambda_amazonas_2)),
  tar_target(df_colombia_sim, create_colombia_sim(name_3, lambda_colombia)),
  tar_target(df_all_muench,
             combine_all_and_pivot(df_amazonas_1, df_amazonas_2, df_colombia,
                                   df_amazonas_sim_1, df_amazonas_sim_2, df_colombia_sim,
                                   name_1, name_2, name_3)),
  tar_target(plot_muench_seroprevalence_val, plot_muench_seroprevalence(df_all_muench)),
  tar_target(plot_muench_foi_val, plot_muench_foi(name_1,
                                              name_2,
                                              name_3,
                                              lambda_amazonas_1,
                                              lambda_amazonas_2,
                                              lambda_colombia)),
  tar_target(plot_muench_both, {
    plot_muench_seroprevalence_val / plot_muench_foi_val +
      plot_annotation(tag_levels = 'A')
  }),
  tar_target(file_plot_muench_both, {
    filename <- "outputs/muench_yf.pdf"
    ggsave(filename, plot_muench_both, width = 9, height = 7)
  })
  
  )

