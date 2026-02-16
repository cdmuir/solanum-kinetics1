# Decompose variance in response variables to phylogenetic, population, and residual components
source("r/header.R")

fit_amphi = read_rds("objects/best_amphi_model.rds")

# Variance decomposition
df_var = fit_amphi |>
  as_draws_df() |>
  select(starts_with("."), starts_with("sd_"), starts_with("sigma")) |>
  rename_with(.fn = \(.x) {
    str_replace(.x, "sigma_", "sd_resid__")
  }, .cols = starts_with("sigma_")) |>
  rename_with(.fn = \(.x) {
    str_remove(.x, "_Intercept")
  }, .cols = ends_with("_Intercept")) |>
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
  split(~ trait) |>
  map(summarize_draws) |>
  map(filter, variable != "trait") |>
  imap_dfr(\(.x, .y) {
    .x |> mutate(trait = .y)
  }) |>
  select(variable, median, q5, q95, trait) |>
  mutate(across(where(is_double), ~ if_else(.x == 0, NA_real_, .x))) |>
  filter(variable != "total_var") |>
  mutate(
    vc = fct_recode(
      variable,
      phylogenetic = "phy",
      population = "accession",
      `within-individual` = "accid",
      `between-individual` = "resid"
    ) |>
      fct_relevel(
        c(
          "phylogenetic",
          "population",
          "between-individual",
          "within-individual"
        )
      ),
    trait1 = fct_recode(
      trait,
      `$\\log(\\lambda)$` = "loglambdamean",
      `$\\log(\\tau)$` = "logtaumean",
      `$\\log(l_\\mathrm{gc})$` = "loggcl",
      `$\\log(f_\\mathrm{gmax})$` = "logfgmax"
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


# Alternate - point range
gp2 =  df_var |>
  ggplot(aes(
    trait1,
    median,
    color = vc,
    ymin = q5,
    ymax = q95
  )) +
  geom_pointrange(position = position_dodge2(width = 0.5)) +
  scale_x_discrete(labels = label_parse()) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_viridis_d() +
  labs(y = "proportion of variance", color = "variance component") +
  theme(axis.title.x = element_blank())

# ggsave("figures/variance.pdf", plot = gp1, width = 6, height = 4)
# ggsave("figures/variance-alt.pdf", plot = gp2, width = 6, height = 4)

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
    `95\\% CI` = glue("[{formatC(q5 * 100, format = 'f', digits = 1)}, {formatC(q95 * 100, format = 'f', digits = 1)}]"),
  ) |>
  filter(`\\% variance` != "NA") |>
  arrange(Trait, vc) |>
  select(Trait, component = vc, `\\% variance`, `95\\% CI`) |>
  write_rds("objects/tbl-variance.rds")
