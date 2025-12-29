# NOTE, I am in process of rejiggering data analysis so I can test model using gs scaled to gmax. This script is different than 07_join-data because that joins the brms summary to stomata anatomy.
# 
# Join stomatal kinetic data with stomatal anatomical data
source("r/header.R")

left_join(
  read_rds("data/rh_curves.rds") |>
    select(acc, id, curve_type, t_sec, gsw, A, curve),
  read_rds("data/stomata.rds") |>
    select(acc = accession, id = replicate, lower_gmax, upper_gmax) |>
    # Calculate total gmax and divide by 1e3 to have same units as gsw
    mutate(
      total_gmax = (lower_gmax + upper_gmax) / 1e3,
      lower_gmax = lower_gmax / 1e3
    ) |>
    select(-upper_gmax) |>
    pivot_longer(
      cols = c(lower_gmax, total_gmax),
      names_to = "gmax_type",
      values_to = "gmax"
    ) |>
    mutate(
      curve_type = case_when(
        gmax_type == "lower_gmax" ~ "1-sided RH",
        gmax_type == "total_gmax" ~ "2-sided RH"
      )
    ),
  by = join_by(acc, id, curve_type)
) |>
  write_rds("data/joined-data.rds")
