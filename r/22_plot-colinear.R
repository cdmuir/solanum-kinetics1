# Figure showing colinearity between phylogenetic and regressione effects of gcl
# on tau

source("r/header.R")

fit = read_rds("objects/fits.rds")$fit[[5]] |>
  as_draws_df() |>
  select(cor_phy__logtaumean_Intercept__loggcl_Intercept,
         b_logtaumean_loggcl)

p = ggplot(
  fit,
  aes(x = cor_phy__logtaumean_Intercept__loggcl_Intercept, y = b_logtaumean_loggcl)
) +
  geom_point(alpha = 0.5) +
  labs(x = "Phylogenetic correlation between $l_\\mathrm{gc}$ on $\\tau$", y = "Regression coffecient of $l_\\mathrm{gc}$ on $\\tau$")

tikz(
  "figures/colinear.tex",
  standAlone = TRUE,
  width = 4,
  height = 4
)
print(p)
dev.off()

system("cd figures; pdflatex colinear.tex; rm colinear.aux colinear.log")
