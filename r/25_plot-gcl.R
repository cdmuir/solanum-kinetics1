# Plot comparing adaxial and abaxial guard cell length among populations
source("r/header.R")

stomata = read_rds("data/stomata.rds")
plant_info = read_rds("data/plant_info.rds") |>
  mutate(light_treatment = fct_recode(light_treatment,
   sun = "high",
   shade = "low"
  )) |>
  select(accession, replicate, light_treatment)

df_stomata = stomata |>
  summarize(across(
    matches("(low|upp)er_guard_cell_length_um"),
    \(.x) median(.x, na.rm = TRUE)
  ), .by = c(accession, replicate)) |>
  left_join(plant_info, by = join_by(accession, replicate)) |>
  summarize(across(
    matches("(low|upp)er_guard_cell_length_um"),
    \(.x) mean(.x, na.rm = TRUE)
  ), .by = c(accession, light_treatment))

ggplot(
  df_stomata,
  aes(
    lower_guard_cell_length_um,
    upper_guard_cell_length_um,
    color = light_treatment
  )
) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  geom_point() +
  scale_x_continuous(limits = c(15, 32.5)) +
  scale_y_continuous(limits = c(15, 32.5)) +
  scale_color_manual(name = "growth\nlight intensity", values = c("sun" = "grey", "shade" = "black")) +
  labs(
    x = expression(Abaxial ~ guard ~ cell ~ length ~ (paste(mu, "m"))),
    y = expression(Adaxial ~ guard ~ cell ~ length ~ (paste(mu, "m")))
  ) +
  coord_equal()

ggsave("figures/gcl.pdf", width = 5, height = 4)
