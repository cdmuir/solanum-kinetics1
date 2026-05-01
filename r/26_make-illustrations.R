# Make stomatal illustrations for figures
# Based on illustrations in Zavala-Paez et al. (2026)
source("r/header.R")

gc_illustrations <- read.ijzip("objects/guard-cell-illustration.zip") |>
  imap_dfr(\(.x, .y) {
    .x$coords |>
    close_polygon() |>
      as.matrix() |>
      list() |>
      st_polygon() |>
      st_sfc() |>
      st_sf() |>
      smooth(method = "ksmooth", smoothness = 3) |>
      st_coordinates() |>
      as_tibble() |>
      mutate(name = .y)
}) |>
  mutate(name = factor(name, levels = c("outer", "inner")),
         x = X - max(X),
         y = Y - max(Y)) |>
  reframe(
    side = rep(c("left", "right"), each = length(x)),
    x = c(x, -x),
    y = c(y, y),
    name = c(name, name)
  ) |>
  rescale_illustration() 

write_rds(gc_illustrations, "objects/gc_illustrations.rds")
