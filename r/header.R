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

convergence_criteria = list(
  rhat_max = 1.05,
  ess_min = 400,
  n_divergent = 10
)

# There was a pretty clear break in the distribution of logtau_mean around 7.
# The four curves above that showed unusual patterns of stomatal closure based
# on visual inspection.

logtau_threshold = 7