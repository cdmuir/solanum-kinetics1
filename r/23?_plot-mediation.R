# Test whether plasticity in stomatal anatomy mediates the effects of light 
# and curve type on tau.

# NOTE: We can't include phylogenetic and residual covariances in mediation 
# because effects of lighttreatment and lightintensity are not varying at those 
# levels.

source("r/header.R") 

# Posterior samples of model parameters
post = read_rds("objects/best_model.rds") |>
  as_draws_df() |>
  select(starts_with("."), starts_with("b_"))

# List of tables of direct and indirect effects for text output
mediation = list(
  # total effect = direct effect of sun on tau + 
  #   mediated effect of sun on tau through fgmax +
  #   mediated effect of sun on tau through gcl
  lighttreatment = post |>
    mutate(
      direct_effect = b_logtaumean_lighttreatmentsun,
      mediated_effect_fgmax = b_logitfgmax_lighttreatmentsun * b_logtaumean_logitfgmax,
      mediated_effect_gcl = b_loggcl_lighttreatmentsun * b_logtaumean_loggcl,
      .keep = "none"
    ),
  
  # total effect = direct effect of high light on tau + 
  #   mediated effect of high light on tau through fgmax 
  # GCL did not change between measurement light intensity treatments
  lightintensity = post |>
    mutate(
      direct_effect = b_logtaumean_lightintensityhigh,
      mediated_effect_fgmax = b_logitfgmax_lightintensityhigh * b_logtaumean_logitfgmax,
      .keep = "none"
    ),
  
  # total effect = direct effect of pseudohypo on tau + 
  #   mediated effect of pseudohypo on tau through fgmax +
  #   mediated effect of pseudohypo on tau through gcl
  curvetype = post |>
    mutate(
      direct_effect = b_logtaumean_curve_typepseudohypo,
      mediated_effect_fgmax = b_logitfgmax_curve_typepseudohypo * b_logtaumean_logitfgmax,
      mediated_effect_gcl = b_loggcl_curve_typepseudohypo * b_logtaumean_loggcl,
      .keep = "none"
    )
  
) |>
  map(\(.x) summarize_draws(.x, estimate = median, quantile2, .args = list(probs = c(0.025, 0.975))))

# Plotting mediation paths as a path diagram
df_edges = post |>
  summarize_draws(estimate = median, quantile2, .args = list(probs = c(0.025, 0.975))) |>
  mutate(
    variable = str_replace(variable, "curve_type", "curvetype")
  ) |>
  separate_wider_delim(variable, delim = "_", names = c("b", "to", "from")) |>
  select(-b) |>
  filter(to %in% c("logtaumean", "logitfgmax", "loggcl"), 
         !(from %in% c("Intercept", "lightintensityhigh", "curvetypepseudohypo"))
  ) |>
  mutate(
    abs_est = abs(estimate),
    sig = as.factor(ifelse(`q2.5` > 0 | `q97.5` < 0, "sig.", "n.s.")),
    sign = ifelse(estimate >= 0, "positive", "negative"),
    label = sprintf("%.2f (%.2f,%.2f)", estimate, `q2.5`, `q97.5`)
  ) |>
  mutate(across(c(from, to), \(.x) recode(.x, 
    logtaumean = "$\\tau$",
    logitfgmax = "$f_\\text{gmax}$",
    loggcl = "$l_\\text{gc}$",
    lighttreatmentsun = "sun treatment"
  )))

df_nodes = tribble(
  ~name, ~x, ~y,
  "sun treatment", 0.1, 0.5,
  "$f_\\text{gmax}$", 0.4, 0.75,
  "$l_\\text{gc}$", 0.4, 0.25,
  "$\\tau$", 0.8, 0.5
)

# join node coordinates to edges for plotting
edge_plot = df_edges |>
  left_join(df_nodes, by = c("from" = "name")) |>
  rename(xstart = x, ystart = y) |>
  left_join(df_nodes, by = c("to" = "name")) |>
  rename(xend = x, yend = y)

# Plot
p <- ggplot() +
  # arrows (use geom_curve for curved mediated paths)
  geom_curve(
    data = edge_plot |> filter(!(from == "sun treatment" & to == "$\\tau$")), # mediated edges curved
    aes(x = xstart, y = ystart, xend = xend, yend = yend,
        color = sign, size = abs_est, linetype = sig),
    curvature = 0.25,
    arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
    lineend = "round",
    show.legend = TRUE
  ) +
  # direct path as a straight line
  geom_segment(
    data = edge_plot |> filter(from == "sun treatment" & to == "$\\tau$"),
    aes(x = xstart, y = ystart, xend = xend, yend = yend,
        color = sign, size = abs_est, linetype = sig),
    arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
    lineend = "round"
  ) +
  # node points
  geom_point(data = df_nodes, aes(x = x, y = y), size = 8, shape = 21, fill = "white") +
  geom_text(data = df_nodes, aes(x = x, y = y, label = name), size = 4, fontface = "bold") +
  # effect labels halfway along each edge
  geom_label_repel(
    data = edge_plot,
    aes(x = (xstart + xend)/2, y = (ystart + yend)/2, label = label),
    size = 3,
    nudge_y = 0.03,
    segment.size = 0
  ) +
  # scales
  scale_size_continuous(range = c(0.5, 2.5), guide = "none") +
  scale_color_manual(values = c(positive = "steelblue", negative = "firebrick")) +
  scale_linetype_manual(values = c(`sig.` = "solid", `n.s.` = "dashed"), guide = guide_legend("significance")) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_fixed()

print(p)
    