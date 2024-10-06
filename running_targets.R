if (!require("pacman")) install.packages("pacman")
p_load(renv, targets, visNetwork, rstan, tidyverse, crew)

# install system-wide cmake (unix terminal):
# /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# (echo; echo 'eval "$(/opt/homebrew/bin/brew shellenv)"') >> /Users/ekamau/.zprofile
# eval "$(/opt/homebrew/bin/brew shellenv)"
# brew install cmake 

getwd()

#** RUN ONCE** - to create a _targets.R script file to configure and define the pipeline:
#use_targets()

# define functions in R scripts in R/ folder

# inspect the pipeline:
tar_manifest(fields = all_of("command"))
tar_visnetwork() # also checks changes within the pipeline

# run the pipeline:
tar_make()
tar_outdated()

tar_read(plot_age_FOI_model)
tar_read(plot_time_FOI_model)
tar_read(plot_ev68_fit)
tar_read(plot_ebola_model_result)
tar_read(plot_hiv_model_fig)

# check warnings:
targets::tar_meta(fields = warnings, complete_only = TRUE)
