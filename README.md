# solanum-kinetics

This repository contains source code associated with the manuscript:

Muir CD, WS Lim. Guard cell size and initial conductance influence stomatal closure kinetics. [*bioRxiv*](https://doi.org/10.64898/2026.05.17.725794).

## License

- **Code** (the `r/` scripts, `Makefile`, and `r/install-packages.R`) is licensed under the [MIT License](LICENSE).
- **Manuscript text, figures, and tables** (the `ms/`, `figures/`, and `tables/` directories) are licensed under [CC BY 4.0](LICENSE-CC-BY).
- **Data** in the `data/` directory are downloaded from [cdmuir/solanum-aa](https://github.com/cdmuir/solanum-aa) and remain subject to that repository's own terms.

## Author contributions

- [Chris Muir](https://cdmuir.netlify.app): Conceptualization, Methodology, Investigation, Visualization, Funding acquisition, Writing
- Wei Shen Lim: Investigation, Writing – review & editing

## Contents

This repository has the following file folders:

- `data`: data files downloaded from [cdmuir/solanum-aa](https://github.com/cdmuir/solanum-aa)
- `figures`: figures generated from _R_ code
- `ms`: manuscript input (e.g. `ms.qmd` and `solanum-kinetics.bib`) and output (`ms.pdf`) files
- `objects`: saved _R_ objects generated from _R_ code
- `r`: _R_ scripts for all data processing and analysis
- `tables`: tables generated from _R_ code

## Prerequisites

To run code and render the manuscript:

- [_R_](https://cran.r-project.org/) version ≥ 4.5.0 and [_RStudio_](https://www.posit.co/) (recommended)
- [Quarto](https://quarto.org/): for rendering `ms/ms.qmd`
- [LaTeX](https://www.latex-project.org/): install the full version or use [**tinytex**](https://yihui.org/tinytex/)
- [GNU Make](https://www.gnu.org/software/make/): type `make --version` in a terminal to check if it is already installed

Before running scripts, install the required _R_ packages:

```r
source("r/install-packages.R")
```

Fitting **brms** models also requires a working CmdStan installation:

```r
cmdstanr::install_cmdstan()
```

## Downloading data and code

1. Download or clone this repository to your machine:

```
git clone git@github.com:cdmuir/solanum-kinetics.git
```

2. Open `solanum-kinetics.Rproj` in [RStudio](https://www.posit.co/)

## Rendering and reproducing results

All steps are coordinated by `make`. Three targets are available:

| Target | What it does |
|--------|-------------|
| `make pdf` | Renders `ms/ms.pdf` from existing computed outputs — no R scripts are run. Use this when outputs are already up to date. |
| `make fast` | Reruns all R scripts **except** the slow brms model fits and other computationally intensive scripts (see **Note**), then renders the manuscript. Requires `objects/weibull/` and `objects/fits/` to already exist from a prior `make all`. |
| `make all` | Runs **every** R script including the brms model fits, then renders. This will take a very long time. |

> **Note:** Nine scripts are slow enough to skip in `make fast` (their pre-computed outputs are shipped instead): `r/02_fit-weibull.R`, `r/03_refit-weibull.R`, and `r/10_fit-all.R` fit Bayesian models with [**brms**](https://paul-buerkner.github.io/brms/) and take many hours; `r/30_refit-vpd.R` is a comparatively quick single brms refit; `r/33_simulate-null.R` runs a 1,000-replicate Monte Carlo simulation; `r/04_calc-r2.R`, `r/05_summarize-pars.R`, `r/06_compare-gsw.R`, and `r/07_plot-curves.R` each iterate over every one of the ~2,100 individual weibull curve fits (over a GB on disk) -- `r/07_plot-curves.R` additionally computes a posterior prediction per curve, making it the slowest of the four -- which is slow I/O/computation rather than a model fit but still too slow for `make fast`.

### Typical workflow

Run the full pipeline once (overnight):

```
make all
```

On subsequent runs, rerun only the fast scripts after making changes:

```
make fast
```

To re-render the manuscript without rerunning any R scripts:

```
make pdf
```

To remove stamp files and the rendered PDF (does not delete data or model objects):

```
make clean
```
