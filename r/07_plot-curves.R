# Plot curve fits
source("r/header.R")

sk_dir1 = "objects/weibull"

accession_info = read_rds("data/accession_info.rds")
plant_info = read_rds("data/plant_info.rds")

r2 = read_rds("objects/r2.rds")

# set plotting parameters
fig_width = 6
fig_height = 5

plan(multisession, workers = 19)

rh_curves = list.files(sk_dir1) |>
  
  future_map(\(.x) {
    id1 = str_remove(.x, "\\.rds$")
    
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
    
    fit_weibull = read_rds(file.path(sk_dir1, .x))
    
    br2 = r2 |>
      filter(id == id1) |>
      mutate(
        s = paste0('Bayes~italic(R)^2=="', sprintf('%.3f', Estimate), '"'),
        t_sec = Inf,
        gsw = Inf
      )
    
    df_new = tibble(t_sec = seq(
      min(fit_weibull$data$t_sec),
      max(fit_weibull$data$t_sec),
      length.out = 1e2
    ))
    
    df_pred = bind_rows(
      posterior_epred(fit_weibull, new = df_new) |>
        as_draws() |>
        summarize_draws() |>
        bind_cols(df_new) |>
        rename(gsw = median)
    )
    
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
    
    ggplot(fit_weibull$data, aes(t_sec, gsw)) +
      geom_line(data = df_pred, color = "darkgrey") +
      geom_point() +
      geom_text(data = br2, mapping = aes(label = s), parse = TRUE, vjust = 1.05,
                hjust = 1.05) +
      labs(
        title = main,
        subtitle = sub,
        x = "time (s)",
        y = expression(italic(g)[sw] ~ (mol ~ m^{
          -2
        } ~ s^{
          -1
        }))
      )
      
    
  }, .progress = TRUE)

pdf("figures/rh-curves.pdf", width = fig_width, height = fig_height)
pb <- progress_bar$new(
  total = length(rh_curves),
  format = "  plotting curves [:bar] :percent eta: :eta",
  clear = FALSE,
  width = 60
)

for (i in seq_along(rh_curves)) {
  print(rh_curves[[i]])
  pb$tick()
}

dev.off()
