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
      # mediated_effect_gcl = b_loggcl_lighttreatmentsun * b_logtaumean_loggcl,
      proportion_mediated = mediated_effect_fgmax / (direct_effect + mediated_effect_fgmax),
      .keep = "none"
    ),
  
  # total effect = direct effect of high light on tau +
  #   mediated effect of high light on tau through fgmax
  # GCL did not change between measurement light intensity treatments
  lightintensity = post |>
    mutate(
      direct_effect = b_logtaumean_lightintensityhigh,
      mediated_effect_fgmax = b_logitfgmax_lightintensityhigh * b_logtaumean_logitfgmax,
      proportion_mediated = mediated_effect_fgmax / (direct_effect + mediated_effect_fgmax),
      .keep = "none"
    ),
  
  # total effect = direct effect of pseudohypo on tau +
  #   mediated effect of pseudohypo on tau through fgmax +
  #   mediated effect of pseudohypo on tau through gcl
  curvetype = post |>
    mutate(
      direct_effect = b_logtaumean_curve_typepseudohypo,
      mediated_effect_fgmax = b_logitfgmax_curve_typepseudohypo * b_logtaumean_logitfgmax,
      # mediated_effect_gcl = b_loggcl_curve_typepseudohypo * b_logtaumean_loggcl,
      proportion_mediated = mediated_effect_fgmax / (direct_effect + mediated_effect_fgmax),
      .keep = "none"
    )
  
) |>
  map(\(.x) summarize_draws(
    .x,
    estimate = median,
    quantile2,
    .args = list(probs = c(0.025, 0.975))
  ))

# Plotting mediation paths as a path diagram

## Edges ----

post_summary = post |>
  summarize_draws(estimate = median, se = \(.x, ...) sd(.x), quantile2, .args = list(probs = c(0.025, 0.975))) |>
  mutate(variable = str_replace(variable, "curve_type", "curvetype")) |>
  separate_wider_delim(variable,
                       delim = "_",
                       names = c("b", "to", "from")) |>
  select(-b) 

### Growth light intensity
df_edges_sun = post_summary |>
  filter(to %in% c("logtaumean", "logitfgmax"), !(
    from %in% c("Intercept", "lightintensityhigh", "curvetypepseudohypo")
  )) |>
  prepare_edges() 

### Measurement light intensity
df_edges_high = post_summary |>
  filter(to %in% c("logtaumean", "logitfgmax"), !(
    from %in% c("Intercept", "lighttreatmentsun", "curvetypepseudohypo")
  )) |>
  prepare_edges()

### Pseudohypo curve type
df_edges_pseudohypo = post_summary |>
  filter(to %in% c("logtaumean", "logitfgmax"), !(
    from %in% c("Intercept", "lightintensityhigh", "lighttreatmentsun")
  )) |>
  prepare_edges()

## Nodes ----
df_nodes = crossing(
  name = c("sun treatment", "high light", "pseudohypo"),
  x = 0.1,
  y = 0.5,
  explanatory = TRUE
) |>
  add_row(name = "$f_\\mathrm{gmax}$", x = 0.4, y = 0.75, explanatory = FALSE) |>
  # add_row(name = "$l_\\mathrm{gc}$", x = 0.4, y = 0.25, explanatory = FALSE) |>
  add_row(name = "$\\tau$", x = 0.8, y = 0.5, explanatory = FALSE)

df_node_labels = df_nodes |>
  crossing(
    treatment = c("sun treatment", "high light", "pseudohypo")
  ) |>
  filter((name == treatment) | !explanatory) |>
  mutate(treatment = case_when(
    treatment == "sun treatment" ~ "Growth\nlight intensity",
    treatment == "high light" ~ "Measurement\nlight intensity",
    treatment == "pseudohypo" ~ "Curve type",
    TRUE ~ treatment
  ))

df_facet_labels = df_node_labels |>
  summarize(treatment = unique(treatment), .by = treatment) |>
  mutate(x = -Inf, y = Inf)

## Join nodes to edges ----

edge_plot = join_nodes_edges(df_edges_sun, df_nodes) |>
  full_join(
    tribble(
      ~from, ~to, ~pstart, ~pend,
      "sun treatment", "$f_\\mathrm{gmax}$", 0.2, 0.8,
      "sun treatment", "$\\tau$", 0.25, 0.9,
      "$f_\\mathrm{gmax}$", "$\\tau$", 0.15, 0.9
    ), by = join_by(to, from)
  ) |>
  mutate(
    xstart1 = xstart + pstart * (xend - xstart),
    ystart1 = ystart + pstart * (yend - ystart),
    xend1 = xstart + pend * (xend - xstart),
    yend1 = ystart + pend * (yend - ystart),
    treatment = "Growth\nlight intensity"
  ) |>
  bind_rows(
    join_nodes_edges(df_edges_high, df_nodes) |>
      full_join(
        tribble(
          ~from, ~to, ~pstart, ~pend,
          "high light", "$f_\\mathrm{gmax}$", 0.2, 0.8,
          "high light", "$\\tau$", 0.2, 0.9,
          "$f_\\mathrm{gmax}$", "$\\tau$", 0.15, 0.9
        ), by = join_by(to, from)
      ) |>
      mutate(
        xstart1 = xstart + pstart * (xend - xstart),
        ystart1 = ystart + pstart * (yend - ystart),
        xend1 = xstart + pend * (xend - xstart),
        yend1 = ystart + pend * (yend - ystart),
        treatment = "Measurement\nlight intensity"
      )
    ) |>
  bind_rows(
    join_nodes_edges(df_edges_pseudohypo, df_nodes) |>
      full_join(
        tribble(
          ~from, ~to, ~pstart, ~pend,
          "pseudohypo", "$f_\\mathrm{gmax}$", 0.2, 0.8,
          "pseudohypo", "$\\tau$", 0.2, 0.9,
          "$f_\\mathrm{gmax}$", "$\\tau$", 0.15, 0.9
        ), by = join_by(to, from)
      ) |>
      mutate(
        xstart1 = xstart + pstart * (xend - xstart),
        ystart1 = ystart + pstart * (yend - ystart),
        xend1 = xstart + pend * (xend - xstart),
        yend1 = ystart + pend * (yend - ystart),
        treatment = "Curve type")
    ) |>
  mutate(angle = 180 / pi * atan2(yend - ystart, xend - xstart))

## Plot ----
p <- ggplot(data = edge_plot) +
  facet_grid(rows = vars(treatment), switch = "y") +
  geom_curve(
    aes(
      x = xstart1,
      y = ystart1,
      xend = xend1,
      yend = yend1,
      color = sign,
      linewidth = effect_size,
      linetype = sig
    ),
    curvature = 0,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed"),
    lineend = "round",
    show.legend = TRUE
  ) +
  geom_label(
    data = df_node_labels,
    mapping = aes(x = x, y = y, label = name),
    fill = "white",
    fontface = 2
  ) +
  geom_text(
    data = df_facet_labels,
    mapping = aes(x = x, y = y, label = treatment),
    fontface = 3,
    hjust = 0,
    vjust = 1
  ) +
  geom_text(
    data = edge_plot,
    aes(
      x = (xstart + xend) / 2,
      y = (ystart + yend) / 2,
      label = label,
      angle = angle
    ),
    size = 3, vjust = -2
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0.4, 0.85)) +
  scale_color_manual(values = c(
    positive = "steelblue",
    negative = "firebrick",
    guide = guide_legend(override.aes = list(linewidth = 100))
  )) +
  scale_linewidth_continuous(
    breaks = c(5, 15, 45),
    name = "effect size") +
  scale_linetype_discrete(guide = if (length(unique(edge_plot$sig)) > 1) {
    guide_legend()
  } else {
    "none"
  }
  ) +
  theme_void() +
  theme(
    legend.key.spacing.y = unit(0.5, "cm"),
    legend.key.width = unit(2, "cm"),
    strip.text = element_blank())
print(p)

tikz(
  "figures/mediation.tex",
  standAlone = TRUE,
  width = 6,
  height = 6
)
print(p)
dev.off()

system("cd figures; pdflatex mediation.tex; rm mediation.aux mediation.log")

write_rds(mediation, "objects/mediation.rds")
