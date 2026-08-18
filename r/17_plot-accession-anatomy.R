# --- Plot anatomical parameters (gcl, fgmax) for each accession by treatment ---

source("r/header.R")

selected_model = read_rds("objects/selected_model.rds")

## Effects of growth light intensity and leaf type on gcl
df1 = selected_model$data |>
  summarize(
    loggcl = mean(loggcl),
    .by = c("accession", "lighttreatment", "leaftype")
  ) |>
  mutate(gcl = exp(loggcl), .keep = "unused")

gp1 = ggplot(df1, aes(leaftype, gcl)) +
  facet_grid(. ~ lighttreatment, scales = "free_y") +
  geom_line(mapping = aes(group = accession), color = "grey") +
  geom_point() +
  labs(x = "leaf type", y = "$l_\\mathrm{gc}$ ($\\si{\\micro\\meter}$)") +
  scale_y_continuous(breaks = seq(15, 25, 5), limits = c(15, 27.5))

## Effects of growth light intensity, measurement light intensity, and leaf type on fgmax
df2 = selected_model$data |>
  mutate(fgmax = exp(loggi - loggmax)) |>
  summarize(
    fgmax = mean(fgmax),
    .by = c("accession", "lighttreatment", "lightintensity", "leaftype")
  )

gp2 = ggplot(df2, aes(lightintensity, fgmax)) +
  facet_grid(leaftype ~ lighttreatment) +
  geom_line(mapping = aes(group = accession), color = "grey") +
  geom_point() +
  scale_y_continuous(breaks = seq(0.05, 0.25, 0.05)) +
  labs(x = "measurement light intensity", y = "$f_\\mathrm{gmax}$ (unitless)")

gp = plot_grid(
  gp1,
  gp2,
  nrow = 2,
  rel_heights = c(0.4, 0.6),
  labels = "auto"
)

options(
  tikzLatexPackages = c(
    getOption("tikzLatexPackages"),
    "\\usepackage{siunitx}"
  )
)

tikz(
  "figures/accession-anatomy.tex",
  standAlone = TRUE,
  width = 6,
  height = 8
)
print(gp)
dev.off()

system("cd figures; pdflatex accession-anatomy.tex; rm accession-anatomy.aux accession-anatomy.log")
