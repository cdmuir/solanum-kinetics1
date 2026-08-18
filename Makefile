# Makefile for solanum-kinetics
#
# Usage:
#   make pdf   -- render ms/ms.pdf from existing computed outputs (fastest)
#   make fast  -- rerun all non-brms R scripts then render; requires brms
#                 outputs (objects/weibull/ and objects/fits/) to already
#                 exist from a prior `make all`
#   make all   -- run every R script (r/00_ through r/38_) including brms
#                 models, then render (slow)
#
# This reflects the public repo layout: r/00_ through r/38_ (the archived
# scripts under r/archive/ -- including the former r/27_plot-conceptual-BD.R,
# which is not part of the manuscript -- are not part of this pipeline),
# plus the data/, figures/, objects/, and tables/ directories those scripts
# read and write.
#
# brms scripts (slow): r/02_fit-weibull.R, r/03_refit-weibull.R,
# r/10_fit-all.R

.PHONY: all fast pdf clean
.DEFAULT_GOAL := pdf

STAMPS := .stamps

$(STAMPS):
	mkdir -p $@

# --------------------------------------------------------------------------
# R script stamp targets, in pipeline order. Each depends on r/header.R,
# r/functions.R, its own script, and the stamps of whichever earlier
# scripts produce the data/objects/tables it reads.
# --------------------------------------------------------------------------

$(STAMPS)/00_load-data: r/00_load-data.R r/header.R r/functions.R | $(STAMPS)
	@echo "==> r/00_load-data.R"
	Rscript r/00_load-data.R
	touch $@

$(STAMPS)/01_join-data: r/01_join-data.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data | $(STAMPS)
	@echo "==> r/01_join-data.R"
	Rscript r/01_join-data.R
	touch $@

# ---- brms: individual Weibull fits (slow) ----
$(STAMPS)/02_fit-weibull: r/02_fit-weibull.R r/header.R r/functions.R \
  $(STAMPS)/01_join-data | $(STAMPS)
	@echo "==> r/02_fit-weibull.R (slow)"
	Rscript r/02_fit-weibull.R
	touch $@

# ---- brms: refit non-converged curves (slow) ----
$(STAMPS)/03_refit-weibull: r/03_refit-weibull.R r/header.R r/functions.R \
  $(STAMPS)/02_fit-weibull | $(STAMPS)
	@echo "==> r/03_refit-weibull.R (slow)"
	Rscript r/03_refit-weibull.R
	touch $@

$(STAMPS)/04_calc-r2: r/04_calc-r2.R r/header.R r/functions.R \
  $(STAMPS)/03_refit-weibull | $(STAMPS)
	@echo "==> r/04_calc-r2.R"
	Rscript r/04_calc-r2.R
	touch $@

$(STAMPS)/05_summarize-pars: r/05_summarize-pars.R r/header.R r/functions.R \
  $(STAMPS)/03_refit-weibull | $(STAMPS)
	@echo "==> r/05_summarize-pars.R"
	Rscript r/05_summarize-pars.R
	touch $@

$(STAMPS)/06_compare-gsw: r/06_compare-gsw.R r/header.R r/functions.R \
  $(STAMPS)/03_refit-weibull | $(STAMPS)
	@echo "==> r/06_compare-gsw.R"
	Rscript r/06_compare-gsw.R
	touch $@

$(STAMPS)/07_plot-curves: r/07_plot-curves.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data $(STAMPS)/03_refit-weibull $(STAMPS)/04_calc-r2 | $(STAMPS)
	@echo "==> r/07_plot-curves.R"
	Rscript r/07_plot-curves.R
	touch $@

$(STAMPS)/08_join-summary: r/08_join-summary.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data $(STAMPS)/01_join-data $(STAMPS)/05_summarize-pars | $(STAMPS)
	@echo "==> r/08_join-summary.R"
	Rscript r/08_join-summary.R
	touch $@

$(STAMPS)/09_make-tbl-vpd: r/09_make-tbl-vpd.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data | $(STAMPS)
	@echo "==> r/09_make-tbl-vpd.R"
	Rscript r/09_make-tbl-vpd.R
	touch $@

# ---- brms: multiresponse phylogenetic model (slow) ----
$(STAMPS)/10_fit-all: r/10_fit-all.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data $(STAMPS)/08_join-summary | $(STAMPS)
	@echo "==> r/10_fit-all.R (slow)"
	Rscript r/10_fit-all.R
	touch $@

$(STAMPS)/11_compare-models: r/11_compare-models.R r/header.R r/functions.R \
  $(STAMPS)/10_fit-all | $(STAMPS)
	@echo "==> r/11_compare-models.R"
	Rscript r/11_compare-models.R
	touch $@

$(STAMPS)/12_select-model: r/12_select-model.R r/header.R r/functions.R \
  $(STAMPS)/11_compare-models | $(STAMPS)
	@echo "==> r/12_select-model.R"
	Rscript r/12_select-model.R
	touch $@

$(STAMPS)/13_plot-estimates: r/13_plot-estimates.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/13_plot-estimates.R"
	Rscript r/13_plot-estimates.R
	touch $@

$(STAMPS)/14_get-partial-cor: r/14_get-partial-cor.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/14_get-partial-cor.R"
	Rscript r/14_get-partial-cor.R
	touch $@

$(STAMPS)/15_make-tbl-estimates-curve: r/15_make-tbl-estimates-curve.R r/header.R r/functions.R \
  $(STAMPS)/09_make-tbl-vpd $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/15_make-tbl-estimates-curve.R"
	Rscript r/15_make-tbl-estimates-curve.R
	touch $@

$(STAMPS)/16_make-tbl-estimates-accession: r/16_make-tbl-estimates-accession.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/16_make-tbl-estimates-accession.R"
	Rscript r/16_make-tbl-estimates-accession.R
	touch $@

$(STAMPS)/17_plot-accession-anatomy: r/17_plot-accession-anatomy.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/17_plot-accession-anatomy.R"
	Rscript r/17_plot-accession-anatomy.R
	touch $@

$(STAMPS)/18_make-tbl-fit-summary: r/18_make-tbl-fit-summary.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/18_make-tbl-fit-summary.R"
	Rscript r/18_make-tbl-fit-summary.R
	touch $@

$(STAMPS)/19_plot-accession-kinetics: r/19_plot-accession-kinetics.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/19_plot-accession-kinetics.R"
	Rscript r/19_plot-accession-kinetics.R
	touch $@

$(STAMPS)/20_plot-variance: r/20_plot-variance.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/20_plot-variance.R"
	Rscript r/20_plot-variance.R
	touch $@

$(STAMPS)/21_plot-gcl-tau: r/21_plot-gcl-tau.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/21_plot-gcl-tau.R"
	Rscript r/21_plot-gcl-tau.R
	touch $@

$(STAMPS)/22_plot-collinear: r/22_plot-collinear.R r/header.R r/functions.R \
  $(STAMPS)/10_fit-all | $(STAMPS)
	@echo "==> r/22_plot-collinear.R"
	Rscript r/22_plot-collinear.R
	touch $@

$(STAMPS)/23_plot-fgmax-kinetics: r/23_plot-fgmax-kinetics.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/23_plot-fgmax-kinetics.R"
	Rscript r/23_plot-fgmax-kinetics.R
	touch $@

$(STAMPS)/24_plot-mediation: r/24_plot-mediation.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/24_plot-mediation.R"
	Rscript r/24_plot-mediation.R
	touch $@

$(STAMPS)/25_plot-gcl: r/25_plot-gcl.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data | $(STAMPS)
	@echo "==> r/25_plot-gcl.R"
	Rscript r/25_plot-gcl.R
	touch $@

$(STAMPS)/26_plot-conceptual: r/26_plot-conceptual.R r/header.R r/functions.R | $(STAMPS)
	@echo "==> r/26_plot-conceptual.R"
	Rscript r/26_plot-conceptual.R
	touch $@

$(STAMPS)/27_summarize-vpd: r/27_summarize-vpd.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data | $(STAMPS)
	@echo "==> r/27_summarize-vpd.R"
	Rscript r/27_summarize-vpd.R
	touch $@

$(STAMPS)/28_plot-vpd-pairs: r/28_plot-vpd-pairs.R r/header.R r/functions.R \
  $(STAMPS)/27_summarize-vpd | $(STAMPS)
	@echo "==> r/28_plot-vpd-pairs.R"
	Rscript r/28_plot-vpd-pairs.R
	touch $@

$(STAMPS)/29_plot-vpd-h2o-rate: r/29_plot-vpd-h2o-rate.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data $(STAMPS)/27_summarize-vpd | $(STAMPS)
	@echo "==> r/29_plot-vpd-h2o-rate.R"
	Rscript r/29_plot-vpd-h2o-rate.R
	touch $@

$(STAMPS)/30_refit-vpd: r/30_refit-vpd.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model $(STAMPS)/27_summarize-vpd | $(STAMPS)
	@echo "==> r/30_refit-vpd.R (slow)"
	Rscript r/30_refit-vpd.R
	touch $@

$(STAMPS)/31_plot-tau-vpd-rate: r/31_plot-tau-vpd-rate.R r/header.R r/functions.R \
  $(STAMPS)/27_summarize-vpd $(STAMPS)/30_refit-vpd | $(STAMPS)
	@echo "==> r/31_plot-tau-vpd-rate.R"
	Rscript r/31_plot-tau-vpd-rate.R
	touch $@

$(STAMPS)/32_plot-comparison: r/32_plot-comparison.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model $(STAMPS)/30_refit-vpd | $(STAMPS)
	@echo "==> r/32_plot-comparison.R"
	Rscript r/32_plot-comparison.R
	touch $@

$(STAMPS)/33_simulate-null: r/33_simulate-null.R r/header.R r/functions.R \
  $(STAMPS)/01_join-data $(STAMPS)/05_summarize-pars | $(STAMPS)
	@echo "==> r/33_simulate-null.R"
	Rscript r/33_simulate-null.R
	touch $@

$(STAMPS)/34_plot-null: r/34_plot-null.R r/header.R r/functions.R \
  $(STAMPS)/08_join-summary $(STAMPS)/33_simulate-null | $(STAMPS)
	@echo "==> r/34_plot-null.R"
	Rscript r/34_plot-null.R
	touch $@

$(STAMPS)/35_validate-nls-vs-bayes: r/35_validate-nls-vs-bayes.R r/header.R r/functions.R \
  $(STAMPS)/01_join-data $(STAMPS)/15_make-tbl-estimates-curve | $(STAMPS)
	@echo "==> r/35_validate-nls-vs-bayes.R"
	Rscript r/35_validate-nls-vs-bayes.R
	touch $@

$(STAMPS)/36_mediation-sensitivity: r/36_mediation-sensitivity.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/36_mediation-sensitivity.R"
	Rscript r/36_mediation-sensitivity.R
	touch $@

$(STAMPS)/37_plot-mediation-dag: r/37_plot-mediation-dag.R r/header.R r/functions.R | $(STAMPS)
	@echo "==> r/37_plot-mediation-dag.R"
	Rscript r/37_plot-mediation-dag.R
	touch $@

$(STAMPS)/38_get-loo-diagnostics: r/38_get-loo-diagnostics.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model | $(STAMPS)
	@echo "==> r/38_get-loo-diagnostics.R"
	Rscript r/38_get-loo-diagnostics.R
	touch $@

# --------------------------------------------------------------------------
# Manuscript
# --------------------------------------------------------------------------

RENDER := cd ms && quarto render ms.qmd

# All stamps, in pipeline order (used by both `fast` and `all`)
NONBRMS_STAMPS := \
  $(STAMPS)/00_load-data \
  $(STAMPS)/01_join-data \
  $(STAMPS)/04_calc-r2 \
  $(STAMPS)/05_summarize-pars \
  $(STAMPS)/06_compare-gsw \
  $(STAMPS)/07_plot-curves \
  $(STAMPS)/08_join-summary \
  $(STAMPS)/09_make-tbl-vpd \
  $(STAMPS)/11_compare-models \
  $(STAMPS)/12_select-model \
  $(STAMPS)/13_plot-estimates \
  $(STAMPS)/14_get-partial-cor \
  $(STAMPS)/15_make-tbl-estimates-curve \
  $(STAMPS)/16_make-tbl-estimates-accession \
  $(STAMPS)/17_plot-accession-anatomy \
  $(STAMPS)/18_make-tbl-fit-summary \
  $(STAMPS)/19_plot-accession-kinetics \
  $(STAMPS)/20_plot-variance \
  $(STAMPS)/21_plot-gcl-tau \
  $(STAMPS)/22_plot-collinear \
  $(STAMPS)/23_plot-fgmax-kinetics \
  $(STAMPS)/24_plot-mediation \
  $(STAMPS)/25_plot-gcl \
  $(STAMPS)/26_plot-conceptual \
  $(STAMPS)/27_summarize-vpd \
  $(STAMPS)/28_plot-vpd-pairs \
  $(STAMPS)/29_plot-vpd-h2o-rate \
  $(STAMPS)/30_refit-vpd \
  $(STAMPS)/31_plot-tau-vpd-rate \
  $(STAMPS)/32_plot-comparison \
  $(STAMPS)/33_simulate-null \
  $(STAMPS)/34_plot-null \
  $(STAMPS)/35_validate-nls-vs-bayes \
  $(STAMPS)/36_mediation-sensitivity \
  $(STAMPS)/37_plot-mediation-dag \
  $(STAMPS)/38_get-loo-diagnostics

BRMS_STAMPS := \
  $(STAMPS)/02_fit-weibull \
  $(STAMPS)/03_refit-weibull \
  $(STAMPS)/10_fit-all

# --------------------------------------------------------------------------
# Top-level targets
# --------------------------------------------------------------------------

# Render ms/ms.pdf using whatever outputs already exist -- no R scripts run
pdf:
	$(RENDER)

# Run all non-brms scripts (requires objects/weibull/ and objects/fits/
# to already exist from a prior `make all`), then render
fast: $(NONBRMS_STAMPS)
	$(RENDER)

# Run every script (r/00_ through r/38_) including brms models, then render
all: $(BRMS_STAMPS) $(NONBRMS_STAMPS)
	$(RENDER)

# --------------------------------------------------------------------------
# Diff: generate a tracked-changes PDF comparing submitted vs revised ms.
# Requires latexdiff (bundled with most TeX distributions).
# Run `make pdf` first so ms/ms.tex is up to date, then run `make diff`.
# --------------------------------------------------------------------------

.PHONY: diff
diff: ms/ms-v1-submitted.tex ms/ms.tex
	latexdiff --encoding=utf8 ms/ms-v1-submitted.tex ms/ms.tex > ms/ms-diff.tex
	cd ms && latexmk -pdf -interaction=nonstopmode ms-diff.tex
	cd ms && latexmk -c ms-diff.tex

# --------------------------------------------------------------------------
clean:
	rm -rf $(STAMPS)
	rm -f ms/ms.pdf
	rm -f ms/ms-diff.tex ms/ms-diff.pdf
