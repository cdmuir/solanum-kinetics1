# Reviewer comment R1.9: document the leaflet-matching fallback in
# r/01_join-data.R and the leaflet-type correction now applied there. The
# correction model itself is fit in r/01_join-data.R (where it is used to
# adjust guard cell length and gmax for individuals whose anatomy leaflet
# did not match their LI-6800 leaflet); this script only reports on it.

source("r/header.R")

# --- Fallback frequency: how often did the LICOR leaflet have no matching
# stomata measurement, requiring a substitute (now leaflet-type-corrected)
# leaflet? ---

plant_info = read_rds("data/plant_info.rds") |>
  filter(!is.na(`1s_rh_response_date`)) |>
  select(acc = accession, id = replicate, leaflet_licor = leaflet)

fallback_check = read_rds("data/stomata.rds") |>
  select(acc = accession, id = replicate, leaflet_stomata = leaflet,
         ends_with("_gmax"), ends_with("_guard_cell_length_um"), ends_with("_stomatal_density_mm2")) |>
  mutate(
    total_gmax = (lower_gmax + upper_gmax) / 1e3,
    lower_gmax = lower_gmax / 1e3,
    total_guard_cell_length_um = (
      lower_guard_cell_length_um * lower_stomatal_density_mm2 + upper_guard_cell_length_um * upper_stomatal_density_mm2
    ) / total_stomatal_density_mm2
  ) |>
  select(-starts_with("upper_"), -ends_with("_stomatal_density_mm2")) |>
  pivot_longer(cols = -c(acc, id, leaflet_stomata), names_to = c("curve_type", "trait"),
               names_pattern = "(lower|upper|total)_(.*)") |>
  mutate(curve_type = case_when(curve_type == "lower" ~ "1-sided RH", curve_type == "total" ~ "2-sided RH")) |>
  pivot_wider(names_from = trait, values_from = value) |>
  full_join(plant_info, by = join_by(acc, id)) |>
  summarize(
    n_match = sum(leaflet_stomata == leaflet_licor & !is.na(gmax), na.rm = TRUE),
    .by = c(acc, id, curve_type)
  )

fallback_summary = list(
  n_total = nrow(fallback_check),
  n_fallback = sum(fallback_check$n_match != 1),
  pct_fallback = 100 * mean(fallback_check$n_match != 1)
)

write_rds(fallback_summary, "objects/gcl-fallback-summary.rds")

# --- Report on the leaflet-type correction model fit in r/01_join-data.R --

# --- Figure: guard cell length by leaflet type, for the curve type (amphi)
# most relevant to the fgmax/gcl analyses. The fitted lm object stores its
# own model frame (leaflet_correction[[ct]]$model$model), so we reuse that
# rather than re-deriving the leaflet-level dataset. ---

amphi_model = leaflet_correction[["2-sided RH"]]$model

p_leaflet_gcl = ggplot(amphi_model$model, aes(leaflet_stomata, exp(`log(guard_cell_length_um)`))) +
  geom_boxplot(outlier.shape = NA, width = 0.4) +
  geom_jitter(width = 0.08, alpha = 0.4) +
  labs(x = NULL, y = expression("Guard cell length (" * mu * "m)"))

ggsave("figures/gcl-leaflet-type.pdf", p_leaflet_gcl, width = 5, height = 4.5)
