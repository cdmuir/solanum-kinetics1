# Plot curve fits
source("r/header.R")

accession_info = read_rds("https://github.com/cdmuir/solanum-aa/raw/refs/heads/main/data/accession-info.rds")
plant_info = read_rds("https://github.com/cdmuir/solanum-aa/raw/refs/heads/main/data/plant-info.rds")

# set plotting parameters
fig_width = 6
fig_height = 5
fontsize_pt = 10
hjust = -0.1
vjust = 0.5


.x = "LA0107-J_pseudohypo_150"

pix = plant_info |>
  filter(acc_id == str_extract(.x, "^(LA[0-9]{4}A*|nelsonii|sandwicense)-[A-Z][AB]{0,1}")) |>
  select(accession, replicate, light_treatment) |>
  mutate(
    curve_type = str_extract(.x, "amphi|pseudohypo"),
    light_intensity = str_extract(.x, "150|2000") |>
      factor(levels = c("150", "2000")) |>
      fct_recode(shade = "150", sun = "2000")
  ) |>
  left_join(select(accession_info, species, accession), by = join_by(accession))

m = read_rds(file.path("objects/sk-curves/", paste0(.x, ".rds")))

br2 = bayes_R2(m)

df_new = tibble(t_sec = seq(min(m$data$t_sec), max(m$data$t_sec), length.out = 1e2))

df_pred = posterior_epred(m, new = df_new) |>
  as_draws() |>
  summarize_draws() |>
  bind_cols(df_new) |>
  rename(gsw = median)

main = glue(
  "{acc}-{rep} ({sp})",
  acc = pix$accession,
  rep = pix$replicate,
  sp = pix$species
)

sub = glue("growth light intensity: {x}",
           x = filter(plant_info, acc_id == .x)$light_treatment)

ggplot(m$data, aes(t_sec, gsw)) +
  geom_line(data = df_pred) +
  geom_point()

summary(m)

bayesplot::mcmc_trace(m)
m


