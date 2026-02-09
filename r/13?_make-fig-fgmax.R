# Plot growth light intensity, measurement light intensity, and leaf treatment
# against fgmax
source("r/header.R")

joined_summary = read_rds("data/joined-summary.rds")

# scratch - checking out plots
ggplot(joined_summary, aes(light_intensity, f_gmax)) +
  geom_jitter() +
  facet_wrap(~light_treatment) +
  scale_y_log10()
