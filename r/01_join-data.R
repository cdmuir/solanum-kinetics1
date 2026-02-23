# Join stomatal kinetic data with stomatal anatomical data
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
  filter(!(curve %in% rh_over5)) |>
  select(acc, id, curve_type, t_sec, gsw, H2O_r, rh_air, curve)

plant_info = read_rds("data/plant_info.rds") |>
  filter(!is.na(`1s_rh_response_date`)) |>
  select(
    acc = accession,
    id = replicate,
    leaflet_licor = leaflet
  )

stomata = read_rds("data/stomata.rds") |>
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
    ) / total_stomatal_density_mm2,
  ) |>
  select(-starts_with("upper_"), -ends_with("_stomatal_density_mm2")) |>
  pivot_longer(
    cols = -c(acc, id, leaflet_stomata),
    names_to = c("curve_type", "trait"),
    names_pattern = "(lower|upper|total)_(.*)"
  ) |>
  mutate(
    curve_type = case_when(
      curve_type == "lower" ~ "1-sided RH",
      curve_type == "total" ~ "2-sided RH"
    )
  ) |>
  pivot_wider(names_from = trait, values_from = value) |>
  full_join(plant_info, by = join_by(acc,  id)) |>
  group_by(acc, id, curve_type) |>
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
      idx_match = which(leaflet_stomata == leaflet_licor & !is.na(guard_cell_length_um))
      if (length(idx_match) == 1) {
        guard_cell_length_um[idx_match[1]] # take the (first) matching value
      } else {
        guard_cell_length_um[which(!is.na(guard_cell_length_um))[1]] # otherwise take first available
      }
    },
    .groups = "drop"
  ) 
     
# Join and write
left_join(
  rh_curves,
  stomata,
  by = join_by(acc, id, curve_type)
) |>
  # Remove those without stomatal anatomy data
  filter(!is.na(gmax)) |>
  # Remove guard cell length outliers (see r/scratch - gcl outliers.R)
  filter(
    !((acc == "LA0407" & id == "P") | 
      (acc == "LA4116" & id == "I") |
      (acc == "LA0429" & id == "F"))
    ) |>
  write_rds("data/joined-data.rds")
