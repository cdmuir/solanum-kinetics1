# quick and dirty model comparison while doing model dev
x = list.files("objects/sk-curves")
i <- 1
m1 = read_rds(paste0("objects/sk-curves/", x[i])) |>
  add_criterion("loo", moment_match = FALSE)
m2 = read_rds(paste0("objects/sk-curves2/", x[i])) |>
  add_criterion("loo", moment_match = FALSE)
m3 = read_rds(paste0("objects/sk-curves3/", x[i])) |>
  add_criterion("loo", moment_match = FALSE)

bayes_R2(m1)
bayes_R2(m2)
bayes_R2(m3)
loo_compare(m1, m2, m3)

conditional_effects(m1)
conditional_effects(m2)

ggplot(m1$data, aes(t_sec, gsw)) +
  geom_point() +
  geom_line(data = conditional_effects(m1)[[1]],
            aes(t_sec, estimate__),
            color = "tomato") +
  geom_ribbon(
    data = conditional_effects(m1)[[1]],
    aes(t_sec, ymin = lower__, ymax = upper__),
    alpha = 0.2,
    fill = "tomato"
  ) +
geom_line(data = conditional_effects(m2)[[1]],
            aes(t_sec, estimate__),
            color = "steelblue") +
  geom_ribbon(
    data = conditional_effects(m2)[[1]],
    aes(t_sec, ymin = lower__, ymax = upper__),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(data = conditional_effects(m3)[[1]],
            aes(t_sec, estimate__),
            color = "red") +
  geom_ribbon(
    data = conditional_effects(m3)[[1]],
    aes(t_sec, ymin = lower__, ymax = upper__),
    alpha = 0.2,
    fill = "red"
  )

