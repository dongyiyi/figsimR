# Basic simulation workflow

## Goal

This tutorial shows how to run a basic fig wasp community simulation.

## Load package

``` r

library(figsimR)
library(igraph)
library(ggraph)
library(patchwork)
library(dplyr)
library(stringr)
library(dplyr)
library(stringr)
library(FSA)        
library(ggplot2)
library(RColorBrewer)
library(purrr)
library(tidyr)
library(ggpubr)
library(vegan)     
library(Hmsc)       # For Hierarchical Modelling of Species Communities (HMSC)
library(reshape2)   
# Load required libraries for parallel computing and progress tracking
library(future)      # for parallel backend
library(progressr)   # to display progress bar during long computation
library(tidyverse)
```

## the default parameter list

``` r

names(parameter_list_default)
parameter_list_default$species_roles
```

## Run a small simulation

``` r

# set seed
set.seed(123)

sim <- simulate_figwasp_community(
  num_figs = 100,
  fecundity_mean = parameter_list_default$fecundity_mean,
  fecundity_dispersion = parameter_list_default$fecundity_dispersion,
  entry_mu = parameter_list_default$entry_mu,
  entry_size = parameter_list_default$entry_size,
  entry_priority = parameter_list_default$entry_priority,
  species_roles = parameter_list_default$species_roles,
  max_entry_table = parameter_list_default$max_entry_table,
  egg_success_prob = parameter_list_default$egg_success_prob,
  egg_success_prob_by_phase = parameter_list_default$egg_success_prob_by_phase,
  layer_preference = parameter_list_default$layer_preference
)

# check the output
names(sim)
head(sim$summary) # the simulated community data
head(sim$original_summary) # more details for wasps and simulations

# Extract emergence matrix (only fig wasp community)

emergence_cols <- grep("^emergence_", names(sim$summary), value = TRUE)
comm_mat <- sim$summary[, emergence_cols]
head(comm_mat)
```
