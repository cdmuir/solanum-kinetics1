source("r/header.R")

curve_vpd = read_rds("objects/curve-vpd-summary.rds")
pairs_data = curve_vpd |>
  select(median_vpd, final_vpd, vpd_rate = sat_rate)

# Pairwise correlations (Pearson)
vpd_cor = cor(pairs_data, use = "pairwise.complete.obs")

# Scatterplot matrix, colored by leaf type since Table S3 highlighted
# leaf-type differences in realized VPD
plot_data = curve_vpd |>
  mutate(
    leaf_type = factor(
      case_when(
        curve_type == "1-sided RH" ~ "pseudohypo",
        curve_type == "2-sided RH" ~ "amphi"
      ),
      levels = c("amphi", "pseudohypo")
    )
  )

p_median_final = ggplot(plot_data, aes(median_vpd, final_vpd, color = leaf_type)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(amphi = col_amphi, pseudohypo = col_pseudohypo)) +
  labs(x = "Median VPD (kPa)", y = "Final VPD (kPa)", color = "Leaf type")

p_median_rate = ggplot(plot_data, aes(median_vpd, sat_rate, color = leaf_type)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(amphi = col_amphi, pseudohypo = col_pseudohypo)) +
  labs(x = "Median VPD (kPa)", y = "VPD saturation rate (1/s)", color = "Leaf type")

p_final_rate = ggplot(plot_data, aes(final_vpd, sat_rate, color = leaf_type)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(amphi = col_amphi, pseudohypo = col_pseudohypo)) +
  labs(x = "Final VPD (kPa)", y = "VPD saturation rate (1/s)", color = "Leaf type")

p_vpd_pairs = (p_median_final + p_median_rate + p_final_rate) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

ggsave("figures/vpd-pairs.pdf", width = 6.5, height = 3)
