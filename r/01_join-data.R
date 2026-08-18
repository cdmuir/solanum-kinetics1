# --- Clean humidity-response curves and join with leaflet-level stomatal anatomy ---
#
# Produces data/joined-data.rds (per-timepoint gsw data joined with gmax, guard
# cell length, and stomatal density) and objects/gcl_switch_summary.rds (how
# often the anatomical leaflet had to be substituted for the LI-6800 leaflet).

source("r/header.R")

rh_curves = read_rds("data/rh_curves.rds")  |>
  mutate(svp = svp(Tair, Pa), rh_air = H2O_r / svp)

rhair = rh_curves |>
  summarize(H2O_r = mean(H2O_r),
            rh_air = mean(rh_air),
            .by = curve)

# Exclude samples where RH_air > 5%
rh_over5 = rhair |>
  filter(rh_air > 5) |>
  pull(curve)

rh_curves = rh_curves |>
  filter(!(curve %in% rh_over5))

# Average time between the start of the high- and low-light-intensity
# curves logged for the same leaf. Curves for the same leaf are typically
# logged within the same file (high light first, then low light after
# re-acclimation); `elapsed` is a continuous absolute time within a file,
# so the gap is the difference in `elapsed` at the start (t_sec == 0) of
# each curve. Exported for use in the Materials and Methods (ms/ms.qmd).

light_intensity_gap = rh_curves |>
  filter(t_sec == 0) |>
  distinct(file, curve, light_intensity, elapsed) |>
  add_count(file) |>
  filter(n == 2) |>
  select(-n) |>
  summarize(start_high = elapsed[light_intensity == "2000"],
            start_low = elapsed[light_intensity == "150"],
            .by = file) |>
  mutate(gap_sec = start_low - start_high) |>
  # filter those that needed to be restarted
  filter(gap_sec > 1000)

mean_light_intensity_gap_min = mean(light_intensity_gap$gap_sec) / 60

write_rds(mean_light_intensity_gap_min,
          "objects/mean-light-intensity-gap-min.rds")

rh_curves = rh_curves |>
  select(acc, id, curve_type, t_sec, gsw, H2O_r, rh_air, curve)

plant_info = read_rds("data/plant_info.rds") |>
  filter(!is.na(`1s_rh_response_date`)) |>
  select(
    acc = accession,
    id = replicate,
    leaflet_licor = leaflet,
    light_treatment
  )

stomata_long = read_rds("data/stomata.rds") |>
  select(
    acc = accession,
    id = replicate,
    leaflet_stomata = leaflet,
    ends_with("_gmax"),
    ends_with("_guard_cell_length_um"),
    ends_with("_stomatal_density_mm2")
  ) |>
  mutate(
    # Calculate total gmax and divide by 1e3 to have same units as gsw
    total_gmax = (lower_gmax + upper_gmax) / 1e3,
    lower_gmax = lower_gmax / 1e3,
    # Calculate weighted average guard cell length
    total_guard_cell_length_um = (
      lower_guard_cell_length_um * lower_stomatal_density_mm2 + upper_guard_cell_length_um * upper_stomatal_density_mm2
    ) / total_stomatal_density_mm2
  ) |>
  select(-starts_with("upper_")) |>
  pivot_longer(
    cols = -c(acc, id, leaflet_stomata),
    names_to = c("curve_type", "trait"),
    names_pattern = "(lower|upper|total)_(.*)"
  ) |>
  mutate(curve_type = case_when(
    curve_type == "lower" ~ "1-sided RH",
    curve_type == "total" ~ "2-sided RH"
  )) |>
  pivot_wider(names_from = trait, values_from = value) |>
  full_join(plant_info, by = join_by(acc, id))

stomata = stomata_long |>
  summarise(
    gmax = {
      idx_match = which(leaflet_stomata == leaflet_licor & !is.na(gmax))
      if (length(idx_match) == 1) {
        gmax[idx_match[1]]                 # take the (first) matching value
      } else {
        gmax[which(!is.na(gmax))[1]]       # otherwise take first available
      }
    },
    guard_cell_length_um = {
      idx_match = which(leaflet_stomata == leaflet_licor &
                          !is.na(guard_cell_length_um))
      if (length(idx_match) == 1) {
        guard_cell_length_um[idx_match[1]] # take the (first) matching value
      } else {
        guard_cell_length_um[which(!is.na(guard_cell_length_um))[1]] # otherwise take first available
      }
    },
    stomatal_density_mm2 = {
      idx_match = which(leaflet_stomata == leaflet_licor & !is.na(stomatal_density_mm2))
      if (length(idx_match) == 1) {
        stomatal_density_mm2[idx_match[1]]                 # take the (first) matching value
      } else {
        stomatal_density_mm2[which(!is.na(stomatal_density_mm2))[1]]       # otherwise take first available
      }
    },
    # Was alternate leaflet used for stomatal anatomy?
    leaflet_switched = {
      idx_match = which(leaflet_stomata == leaflet_licor &
                          !is.na(guard_cell_length_um))
      length(idx_match) != 1
    },
    .by = c(acc, id, curve_type)
  )

joined_data_full = left_join(rh_curves, stomata, by = join_by(acc, id, curve_type)) |>
  # Remove those without stomatal anatomy data
  filter(!is.na(gmax)) |>
  # Remove guard cell length outliers (residual greater than 0.45 in lm(loggcl ~ lighttreatment + leaftype + phy))
  filter(!((acc == "LA0407" & id == "P") |
             (acc == "LA4116" & id == "I") |
             (acc == "LA0429" & id == "F")
  ))

# Number of individual-by-curve-type combinations in the final dataset whose
# guard cell length/gmax were actually leaflet-type-switched (i.e., the
# leaflet substituted for stomatal anatomy did not match the LI-6800
# leaflet).
gcl_switch_summary = list(
  n_switched = joined_data_full |>
    distinct(acc, id, leaflet_switched) |>
    filter(leaflet_switched) |>
    nrow(),
  n_total = joined_data_full |> distinct(acc, id) |> nrow()
)

write_rds(gcl_switch_summary, "objects/gcl_switch_summary.rds")

joined_data_full |>
  select(-leaflet_switched) |>
  write_rds("data/joined-data.rds")
