# Plot kinetic parameters for each accession in each treatment
source("r/header.R")

dat = read_rds("objects/best_model.rds")$data |>
  summarize(
    logtau = mean(logtaumean),
    loglambda = mean(loglambdamean),
    .by = c("accession", "lighttreatment", "lightintensity")) |>
  mutate(
    tau = exp(logtau),
    lambda = exp(loglambda),
    .keep = "unused"
  ) |>
  pivot_longer(
    c(tau, lambda),
    names_to = "parameter",
    values_to = "value"
  ) |>
  mutate(par1 = case_when(
    parameter == "tau" ~ "tau~(s)",
    parameter == "lambda" ~ "lambda~(unitless)"
  ))


ggplot(dat, aes(lightintensity, value)) +
  facet_grid(par1 ~ lighttreatment, scales = "free_y",
             labeller = "label_parsed") +
  geom_line(mapping = aes(group = accession), color = "grey") + 
  geom_point() +
  labs(
    x = "measurement light intensity",
    y = "kinetic trait value (log-scale)") +
  facetted_pos_scales(
    y = list(
      par1 == "tau~(s)" ~ scale_y_log10(breaks = c(100, 200, 400), limits = c(100, 400)),
      par1 == "lambda~(unitless)" ~ scale_y_log10(breaks = seq(1, 1.75, by = 0.25), limits = c(0.95, 1.75))
    )
  )

ggsave("figures/accession-kinetics.pdf", width = 6, height = 6)
