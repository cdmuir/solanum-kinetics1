# Plot anatomical parameters for each accession in each treatment
source("r/header.R")

fit = read_rds("objects/best_model.rds")

## Effects of growth light intensity and curve type on gcl
df1 = fit$data |>
  summarize(
    loggcl = mean(loggcl),
    .by = c("accession", "lighttreatment", "curve_type")) |>
  mutate(
    gcl = exp(loggcl),
    .keep = "unused"
  ) 

gp1 = ggplot(df1, aes(curve_type, gcl)) +
  facet_grid(. ~ lighttreatment, scales = "free_y") +
  geom_line(mapping = aes(group = accession), color = "grey") + 
  geom_point() +
  labs(
    x = expression(curve~type),
    y = expression(paste(guard~cell~length, ' (', mu, 'm)'))
    ) +
  scale_y_continuous(breaks = seq(15, 25, 5), limits = c(15, 27.5))

## Effects of growth light intensity, measurement light intensity, and curve type on fgmax
df2 = fit$data |>
  summarize(
    logitfgmax = mean(logitfgmax),
    .by = c("accession", "lighttreatment", "lightintensity", "curve_type")) |>
  mutate(
    fgmax = plogis(logitfgmax),
    .keep = "unused"
  ) 

gp2 = ggplot(df2, aes(lightintensity, fgmax)) +
  facet_grid(curve_type ~ lighttreatment, scales = "free_y",
             labeller = "label_parsed") +
  geom_line(mapping = aes(group = accession), color = "grey") + 
  geom_point() +
  labs(
    x = "measurement light intensity",
    y = expression(italic(f)[gmax]~(unitless))) #+
  # facetted_pos_scales(
  #   y = list(
  #     par1 == "guard~cell~length~(paste(mu, 's)')" ~ scale_y_log10(breaks = c(100, 200, 400)),
  #     par1 == "italic(f)[gmax]~(unitless)" ~ scale_y_log10(breaks = seq(1, 1.75, by = 0.25), limits = c(0.95, 1.75))
  #   )
  # )

plot_grid(gp1, gp2, nrow = 2, rel_heights = c(0.4, 0.6), align = "hv",
          labels = "auto")
ggsave("figures/accession-anatomy.pdf", width = 6, height = 8)
