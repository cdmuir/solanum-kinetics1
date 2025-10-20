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
library(ggplot2)
library(ggpp)
library(glue)
library(magrittr)
library(posterior)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

theme_set(theme_cowplot())

form_vico = gsw ~ gf + dg * exp(-(t_sec / tau))  
form_cdweibull = gsw ~ gf + dg * exp(-(t_sec / tau) ^ lambda)

convergence_criteria = list(
  rhat_max = 1.05,
  ess_min = 400,
  n_divergent = 10
)
