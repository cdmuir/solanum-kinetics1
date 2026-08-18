# --- Estimate phylogenetic partial correlations among traits -----------------

source("r/header.R")

selected_model = read_rds("objects/selected_model.rds")

draws_df = selected_model |>
  as_draws_df() |>
  select(starts_with("."),
         starts_with("sd_phy__"),
         starts_with("cor_phy__"))

Omega_draws = make_precision_phy(draws_df)

P_post = map_dfr(seq_len(dim(Omega_draws)[1]), \(.i) {
  P = -cov2cor(Omega_draws[.i, , ])
  P1 = P |>
    c() |>
    as.list() |>
    as.data.frame()
  
  names(P1) = outer(rownames(P), colnames(P), paste, sep = "__") |> as.vector()
  
  P1
})

P_post |>
  mutate(.draw = row_number()) |>
  pivot_longer(-.draw) |>
  summarize(
    median = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    .by = name
  ) |>
  write_rds("objects/partial-cor.rds")
