# --- Decompose trait variance into phylogenetic, population, and individual components ---

source("r/header.R")

selected_model = read_rds("objects/selected_model.rds")

# Variance decomposition
df_var = selected_model |>
  as_draws_df() |>
  select(starts_with("."), starts_with("sd_"), starts_with("sigma")) |>
  clean_posterior_names() |>
  pivot_longer(
    cols = -starts_with("."),
    names_sep = "__",
    names_to = c("component", "trait"),
    values_to = "sd"
  ) |>
  mutate(component = str_remove(component, "sd_")) |>
  mutate(var = sd^2, .keep = "unused") |>
  pivot_wider(values_from = var, names_from = component) |>
  mutate(across(where(is_double), ~ if_else(is.na(.x), 0, .x))) |>
  rowwise() |>
  mutate(total_var = sum(c_across(where(is_double)))) |>
  ungroup() |>
  mutate(across(where(is_double), ~ .x / total_var)) |>
  mutate(ind = resid + accid, .keep = "unused") |>
  split(~ trait) |>
  map(summarize_draws, median, quantile2, .args = list(probs = c(0.025, 0.975))) |>
  map(filter, variable != "trait") |>
  imap_dfr(\(.x, .y) {
    .x |> mutate(trait = .y)
  }) |>
  select(variable, median, `q2.5`, `q97.5`, trait) |>
  mutate(across(where(is_double), ~ if_else(.x == 0, NA_real_, .x))) |>
  filter(variable != "total_var") |>
  mutate(
    vc = factor(
      recode_variance_component(variable),
      levels = c("phylogenetic", "population (nonphylogenetic)", "among-individual")
    ),
    trait1 = factor(
      trait_latex_label(trait),
      levels = trait_latex_label(c("loggcl", "loggi", "loggmax", "loglambdamean", "logtaumean"))
    )
  )

# Bar plot
gp1 = df_var |>
  ggplot(aes(median, trait1, fill = vc)) +
  geom_bar(stat = "identity", position = "stack") +
  # scale_y_discrete(labels = label_parse()) +
  scale_fill_viridis_d() +
  labs(x = "proportion of variance", fill = "variance component") +
  theme(axis.title.y = element_blank())

tikz(
  "figures/variance.tex",
  standAlone = TRUE,
  width = 6,
  height = 2
)
print(gp1)
dev.off()

system("cd figures; pdflatex variance.tex; rm variance.aux variance.log")

# Table
df_var |>
  mutate(
    Trait = trait1,
    `\\% variance` = formatC(median * 100, format = "f", digits = 1),
    `95\\% CI` = glue("[{formatC(`q2.5` * 100, format = 'f', digits = 1)}, {formatC(`q97.5` * 100, format = 'f', digits = 1)}]"),
  ) |>
  filter(`\\% variance` != "NA") |>
  arrange(Trait, vc) |>
  select(Trait, component = vc, `\\% variance`, `95\\% CI`) |>
  write_rds("objects/tbl-variance.rds")
