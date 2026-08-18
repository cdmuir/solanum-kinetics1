# --- Plot fixed-effect and phylogenetic-correlation estimates across all plausible models ---

source("r/header.R")

selected_model = read_rds("objects/tbl-comparison.rds") |>
  filter(selected) |>
  pull(model)

estimates_by_model = read_rds("objects/estimates_by_model.rds") |>
  mutate(
    explanatory1 = case_when(
      str_detect(explanatory, "loggcl") ~ trait_latex_label("loggcl"),
      str_detect(explanatory, "loggi") ~ trait_latex_label("loggi"),
      str_detect(explanatory, "loggmax") ~ trait_latex_label("loggmax")
    ),
    response1 = case_when(
      str_detect(response, "logtaumean") ~ trait_latex_label("logtaumean"),
      str_detect(response, "loglambdamean") ~ trait_latex_label("loglambdamean")
    ),
    `significant?` = case_when(
      overlaps_zero | is.na(overlaps_zero) ~ "no/NA",
      !overlaps_zero ~ "yes"
    ),
    s = ifelse(model == selected_model, "*", ""),
    model1 = factor(
      str_c(str_replace(model, "s/fit_", " "), 
            ifelse(model == selected_model, "*", "")), 
      levels = rev(str_c(str_replace(levels(model), "s/fit_", " "), 
      ifelse(levels(model) == selected_model, "*", "")))
    )
  ) |>
  replace_na(list(estimate = 0, q2.5 = 0, q97.5 = 0))

gp1 = ggplot(
  filter(estimates_by_model, type == "b"),
  aes(
    estimate,
    model1,
    xmin = `q2.5`,
    xmax = `q97.5`,
    color = `significant?`
  )
) +
  facet_grid(explanatory1 ~ response1, scales = "free_x") +
  geom_hline(aes(yintercept = model1), color = "grey70", linetype = "dotted") +
  geom_vline(xintercept = 0,
             color = "grey",
             linetype = "dashed") +
  geom_pointinterval(linewidth = 2) +
  scale_color_manual(values = c("grey", "black")) +
  labs(
    x = "fixed effect ($\\pm 95$\\% CIs)"
  ) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "bottom"
  )

gp2 = gp1 %+% filter(estimates_by_model, type == "cor_phy") +
  labs(x = "phylogenetic correlation ($\\pm 95$\\% CIs)") +
  xlim(-1, 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

tikz(
  "figures/estimates-b.tex",
  standAlone = TRUE,
  width = 4,
  height = 6
)
print(gp1)
dev.off()

tikz(
  "figures/estimates-phy.tex",
  standAlone = TRUE,
  width = 4,
  height = 6
)
print(gp2)
dev.off()

system("cd figures; pdflatex estimates-b.tex; pdflatex estimates-phy.tex; rm estimates-b.aux estimates-b.log estimates-phy.aux estimates-phy.log")
