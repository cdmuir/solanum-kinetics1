rm(list = ls())

source("r/functions.R")

library(ape)
library(brms)
library(checkmate)
library(cowplot)
library(dplyr)
library(forcats)
library(furrr)
library(future)
library(ggh4x)
library(ggplot2)
library(ggpp)
library(ggpubr)
library(glue)
library(magrittr)
library(posterior)
library(progress)
library(purrr)
library(readr)
library(scales)
library(stringr)
library(tibble)
library(tidyr)
library(tikzDevice)
library(zip)

theme_set(theme_cowplot())

form_cdweibull = gsw ~ gf + dg * exp(-(t_sec / tau) ^ lambda)

# start convergence_criteria
convergence_criteria = list(
  rhat_max = 1.05,
  ess_min = 400,
  n_divergent = 10
)
# end convergence_criteria
# 
# There was a pretty clear break in the distribution of logtau_mean around 7.
# Four curves above showed unusual patterns of stomatal closure based on visual
# inspection.

logtau_threshold = 7