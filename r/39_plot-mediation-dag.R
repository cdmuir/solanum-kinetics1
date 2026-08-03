# Reviewer comment R1.6: draw a simple DAG illustrating the assumed causal
# structure behind the fgmax path analysis (r/22_plot-mediation.R) and the
# sequential-ignorability assumption tested by the sensitivity analysis in
# r/38_mediation-sensitivity.R: Treatment -> fgmax -> tau, Treatment -> tau
# (direct path), and an unmeasured confounder U with arrows into both fgmax
# and tau (representing, e.g., hydraulic status, measurement order, or
# realized VPD trajectory) -- exactly what the sensitivity parameter rho
# parameterizes.
#
# Built manually with ggplot2 (geom_segment/geom_label), matching the
# hand-rolled path-diagram style already used in r/22_plot-mediation.R,
# rather than adding a new dependency (ggdag/dagitty are not installed).

source("r/header.R")

nodes = tibble(
  label = c("Treatment", "italic(f[gmax])", "italic(tau)", "italic(U)"),
  x = c(0, 0.5, 1, 0.75),
  y = c(0.15, 0.55, 0.15, 1),
  unmeasured = c(FALSE, FALSE, FALSE, TRUE)
)

edges_raw = tibble(
  x = c(0, 0.5, 0, 0.75, 0.75),
  y = c(0.15, 0.55, 0.15, 1, 1),
  xend = c(0.5, 1, 1, 0.5, 1),
  yend = c(0.55, 0.15, 0.15, 0.55, 0.15),
  # dashed edges indicate paths involving the unmeasured confounder U
  dashed = c(FALSE, FALSE, FALSE, TRUE, TRUE)
)

# Shrink each segment along its own direction so arrowheads stop at node
# borders instead of overlapping the node labels.
shrink_segment = function(x, y, xend, yend, inset) {
  dx = xend - x
  dy = yend - y
  len = sqrt(dx ^ 2 + dy ^ 2)
  ux = dx / len
  uy = dy / len
  tibble(
    x = x + ux * inset,
    y = y + uy * inset,
    xend = xend - ux * inset,
    yend = yend - uy * inset
  )
}

edges = map_dfr(seq_len(nrow(edges_raw)), function(i) {
  shrink_segment(
    edges_raw$x[i], edges_raw$y[i], edges_raw$xend[i], edges_raw$yend[i],
    inset = if (edges_raw$dashed[i]) 0.09 else 0.075
  ) |>
    mutate(dashed = edges_raw$dashed[i])
})

p_dag = ggplot() +
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend, linetype = dashed),
    arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
    linewidth = 0.6,
    color = "grey20",
    show.legend = FALSE
  ) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed")) +
  geom_label(
    data = nodes,
    aes(x = x, y = y, label = label, color = unmeasured),
    parse = TRUE,
    size = 5,
    fontface = 2,
    label.padding = unit(0.4, "lines")
  ) +
  scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "grey50"), guide = "none") +
  scale_x_continuous(limits = c(-0.15, 1.15)) +
  scale_y_continuous(limits = c(0, 1.15)) +
  theme_void()

ggsave("figures/mediation-dag.pdf", p_dag, width = 6, height = 4.5)
