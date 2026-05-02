# Makefile for solanum-kinetics1
#
# Usage:
#   make pdf   -- render ms/ms.pdf from existing computed outputs (fastest)
#   make fast  -- rerun all non-brms R scripts then render; requires brms
#                 outputs (objects/weibull/ and objects/fits.rds) to already
#                 exist from a prior `make all`
#   make all   -- run every R script including brms models then render (slow)
#
# brms scripts (slow): r/02_fit-weibull.R, r/03_refit-weibull.R, r/10_fit-all.R

.PHONY: all fast pdf clean
.DEFAULT_GOAL := pdf

STAMPS := .stamps

$(STAMPS):
	mkdir -p $@

# --------------------------------------------------------------------------
# R script stamp targets
# --------------------------------------------------------------------------

$(STAMPS)/00_load-data: r/00_load-data.R r/header.R | $(STAMPS)
	Rscript r/00_load-data.R
	touch $@

$(STAMPS)/01_join-data: r/01_join-data.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data | $(STAMPS)
	Rscript r/01_join-data.R
	touch $@

# ---- brms: individual Weibull fits (slow) ----
$(STAMPS)/02_fit-weibull: r/02_fit-weibull.R r/header.R r/functions.R \
  $(STAMPS)/01_join-data | $(STAMPS)
	Rscript r/02_fit-weibull.R
	touch $@

# ---- brms: refit non-converged curves (slow) ----
$(STAMPS)/03_refit-weibull: r/03_refit-weibull.R r/header.R r/functions.R \
  $(STAMPS)/02_fit-weibull | $(STAMPS)
	Rscript r/03_refit-weibull.R
	touch $@

$(STAMPS)/04_calc-r2: r/04_calc-r2.R r/header.R \
  $(STAMPS)/03_refit-weibull | $(STAMPS)
	Rscript r/04_calc-r2.R
	touch $@

$(STAMPS)/05_summarize-pars: r/05_summarize-pars.R r/header.R \
  $(STAMPS)/03_refit-weibull | $(STAMPS)
	Rscript r/05_summarize-pars.R
	touch $@

$(STAMPS)/06_compare-gsw: r/06_compare-gsw.R r/header.R \
  $(STAMPS)/03_refit-weibull | $(STAMPS)
	Rscript r/06_compare-gsw.R
	touch $@

$(STAMPS)/07_plot-curves: r/07_plot-curves.R r/header.R \
  $(STAMPS)/04_calc-r2 | $(STAMPS)
	Rscript r/07_plot-curves.R
	touch $@

$(STAMPS)/08_join-summary: r/08_join-summary.R r/header.R \
  $(STAMPS)/05_summarize-pars $(STAMPS)/00_load-data | $(STAMPS)
	Rscript r/08_join-summary.R
	touch $@

$(STAMPS)/09_make-tbl-vpd: r/09_make-tbl-vpd.R r/header.R \
  $(STAMPS)/00_load-data | $(STAMPS)
	Rscript r/09_make-tbl-vpd.R
	touch $@

# ---- brms: multiresponse phylogenetic model (slow) ----
$(STAMPS)/10_fit-all: r/10_fit-all.R r/header.R r/functions.R \
  $(STAMPS)/08_join-summary | $(STAMPS)
	Rscript r/10_fit-all.R
	touch $@

$(STAMPS)/11_compare-models: r/11_compare-models.R r/header.R \
  $(STAMPS)/10_fit-all | $(STAMPS)
	Rscript r/11_compare-models.R
	touch $@

$(STAMPS)/12_get-partial-cor: r/12_get-partial-cor.R r/header.R \
  $(STAMPS)/11_compare-models | $(STAMPS)
	Rscript r/12_get-partial-cor.R
	touch $@

$(STAMPS)/13_make-tbl-estimates-curve: r/13_make-tbl-estimates-curve.R r/header.R \
  $(STAMPS)/11_compare-models | $(STAMPS)
	Rscript r/13_make-tbl-estimates-curve.R
	touch $@

$(STAMPS)/14_make-tbl-estimates-accession: r/14_make-tbl-estimates-accession.R r/header.R \
  $(STAMPS)/11_compare-models | $(STAMPS)
	Rscript r/14_make-tbl-estimates-accession.R
	touch $@

$(STAMPS)/15_plot-accession-anatomy: r/15_plot-accession-anatomy.R r/header.R \
  $(STAMPS)/11_compare-models | $(STAMPS)
	Rscript r/15_plot-accession-anatomy.R
	touch $@

$(STAMPS)/16_make-tbl-fit-summary: r/16_make-tbl-fit-summary.R r/header.R \
  $(STAMPS)/11_compare-models | $(STAMPS)
	Rscript r/16_make-tbl-fit-summary.R
	touch $@

$(STAMPS)/17_plot-accession-kinetics: r/17_plot-accession-kinetics.R r/header.R \
  $(STAMPS)/11_compare-models | $(STAMPS)
	Rscript r/17_plot-accession-kinetics.R
	touch $@

$(STAMPS)/18_plot-variance: r/18_plot-variance.R r/header.R \
  $(STAMPS)/11_compare-models | $(STAMPS)
	Rscript r/18_plot-variance.R
	touch $@

$(STAMPS)/19_plot-gcl-tau: r/19_plot-gcl-tau.R r/header.R \
  $(STAMPS)/11_compare-models | $(STAMPS)
	Rscript r/19_plot-gcl-tau.R
	touch $@

$(STAMPS)/20_plot-colinear: r/20_plot-colinear.R r/header.R \
  $(STAMPS)/10_fit-all | $(STAMPS)
	Rscript r/20_plot-colinear.R
	touch $@

$(STAMPS)/21_plot-fgmax-kinetics: r/21_plot-fgmax-kinetics.R r/header.R \
  $(STAMPS)/11_compare-models | $(STAMPS)
	Rscript r/21_plot-fgmax-kinetics.R
	touch $@

$(STAMPS)/22_plot-mediation: r/22_plot-mediation.R r/header.R \
  $(STAMPS)/11_compare-models | $(STAMPS)
	Rscript r/22_plot-mediation.R
	touch $@

$(STAMPS)/23_plot-gcl: r/23_plot-gcl.R r/header.R \
  $(STAMPS)/00_load-data | $(STAMPS)
	Rscript r/23_plot-gcl.R
	touch $@

$(STAMPS)/24_make-illustrations: r/24_make-illustrations.R r/header.R | $(STAMPS)
	Rscript r/24_make-illustrations.R
	touch $@

$(STAMPS)/25_plot-conceptual: r/25_plot-conceptual.R r/header.R \
  $(STAMPS)/24_make-illustrations | $(STAMPS)
	Rscript r/25_plot-conceptual.R
	touch $@

# --------------------------------------------------------------------------
# Manuscript
# --------------------------------------------------------------------------

RENDER := cd ms && quarto render ms.qmd

# All stamps needed before rendering (brms-dependent scripts included)
RENDER_DEPS := \
  $(STAMPS)/04_calc-r2 \
  $(STAMPS)/05_summarize-pars \
  $(STAMPS)/06_compare-gsw \
  $(STAMPS)/07_plot-curves \
  $(STAMPS)/08_join-summary \
  $(STAMPS)/09_make-tbl-vpd \
  $(STAMPS)/11_compare-models \
  $(STAMPS)/12_get-partial-cor \
  $(STAMPS)/13_make-tbl-estimates-curve \
  $(STAMPS)/14_make-tbl-estimates-accession \
  $(STAMPS)/15_plot-accession-anatomy \
  $(STAMPS)/16_make-tbl-fit-summary \
  $(STAMPS)/17_plot-accession-kinetics \
  $(STAMPS)/18_plot-variance \
  $(STAMPS)/19_plot-gcl-tau \
  $(STAMPS)/20_plot-colinear \
  $(STAMPS)/21_plot-fgmax-kinetics \
  $(STAMPS)/22_plot-mediation \
  $(STAMPS)/23_plot-gcl \
  $(STAMPS)/24_make-illustrations \
  $(STAMPS)/25_plot-conceptual

# --------------------------------------------------------------------------
# Top-level targets
# --------------------------------------------------------------------------

# Render ms/ms.pdf using whatever outputs already exist — no R scripts run
pdf:
	$(RENDER)

# Run all non-brms scripts (requires objects/weibull/ and objects/fits.rds
# to already exist from a prior `make all`), then render
fast: \
  $(STAMPS)/00_load-data \
  $(STAMPS)/01_join-data \
  $(STAMPS)/04_calc-r2 \
  $(STAMPS)/05_summarize-pars \
  $(STAMPS)/06_compare-gsw \
  $(STAMPS)/07_plot-curves \
  $(STAMPS)/08_join-summary \
  $(STAMPS)/09_make-tbl-vpd \
  $(STAMPS)/11_compare-models \
  $(STAMPS)/12_get-partial-cor \
  $(STAMPS)/13_make-tbl-estimates-curve \
  $(STAMPS)/14_make-tbl-estimates-accession \
  $(STAMPS)/15_plot-accession-anatomy \
  $(STAMPS)/16_make-tbl-fit-summary \
  $(STAMPS)/17_plot-accession-kinetics \
  $(STAMPS)/18_plot-variance \
  $(STAMPS)/19_plot-gcl-tau \
  $(STAMPS)/20_plot-colinear \
  $(STAMPS)/21_plot-fgmax-kinetics \
  $(STAMPS)/22_plot-mediation \
  $(STAMPS)/23_plot-gcl \
  $(STAMPS)/24_make-illustrations \
  $(STAMPS)/25_plot-conceptual
	$(RENDER)

# Run every script including brms models, then render
all: \
  $(STAMPS)/02_fit-weibull \
  $(STAMPS)/03_refit-weibull \
  $(STAMPS)/10_fit-all \
  $(RENDER_DEPS)
	$(RENDER)

# --------------------------------------------------------------------------
clean:
	rm -rf $(STAMPS)
	rm -f ms/ms.pdf
