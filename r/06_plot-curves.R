# Plot curve fits
source("r/header.R")

accession_info = read_rds("data/accession_info.rds")
plant_info = read_rds("data/plant_info.rds")

r2 = read_rds("objects/r2.rds")

# set plotting parameters
fig_width = 6
fig_height = 5

plan(multisession, workers = 19)

rh_curves = list.files("objects/sk-curves") |>
  str_remove(".rds$") |>
  future_map(\(.x) {
    pix = plant_info |>
      filter(
        acc_id == str_extract(.x, "^(LA[0-9]{4}A*|nelsonii|sandwicense)-[A-Z][AB]{0,1}")
      ) |>
      select(accession, replicate, light_treatment) |>
      mutate(
        curve_type = str_extract(.x, "amphi|pseudohypo"),
        light_treatment = light_treatment |>
          factor(levels = c("low", "high")) |>
          fct_recode(shade = "low", sun = "high"),
        light_intensity = str_extract(.x, "150|2000") |>
          factor(levels = c("150", "2000")) |>
          fct_recode(low = "150", high = "2000")
      ) |>
      left_join(select(accession_info, species, accession), by = join_by(accession))
    
    m = read_rds(file.path("objects/sk-curves/", paste0(.x, ".rds")))
    
    br2 = r2 |>
      filter(id == .x) |>
      pull(Estimate)
    
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
    
    sub = glue(
      "growth: {x}; meas: {y}",
      x = pix$light_treatment,
      y = pix$light_intensity
    )
    
    ggplot(m$data, aes(t_sec, gsw)) +
      geom_line(data = df_pred) +
      geom_point() +
      labs(title = main,
           subtitle = sub,
           x = "time (s)") +
      ylab(expression(italic(g)[sw] ~ (mol ~ m^{
        -2
      } ~ s^{
        -1
      }))) +
      annotate(
        "text_npc",
        npcx = 0.95,
        npcy = 0.95,
        label = paste0('Bayes~italic(R)^2==\"', sprintf('%.3f', br2), '\"'),
        parse = TRUE
      )
  }, .progress = TRUE)

pdf("figures/rh-curves.pdf", width = fig_width, height = fig_height)
for (i in seq_along(rh_curves)) {
  print(rh_curves[[i]])
}

dev.off()

