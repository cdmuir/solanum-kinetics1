# Makefile for solanum-kinetics
#
# Usage:
#   make pdf   -- render ms/ms.pdf from existing computed outputs (fastest)
#   make fast  -- rerun all fast R scripts then render; requires the slow
#                 scripts' outputs (objects/weibull/, objects/fits/,
#                 objects/r2.rds, objects/pars-summary.rds,
#                 figures/compare-gsw.pdf, figures/rh-curves.pdf,
#                 objects/selected_model_vpd.rds, and
#                 objects/null-sim-fgmax-tau.rds) to already exist
#                 from a prior `make all`
#   make all   -- run every R script (r/00_ through r/38_) including the
#                 slow ones, then render (slow)
#
# This reflects the public repo layout: r/00_ through r/38_ (the archived
# scripts under r/archive/ -- including the former r/27_plot-conceptual-BD.R,
# which is not part of the manuscript -- are not part of this pipeline),
# plus the data/, figures/, objects/, and tables/ directories those scripts
# read and write.
#
# Slow scripts skipped by `make fast`: r/02_fit-weibull.R,
# r/03_refit-weibull.R, r/10_fit-all.R, and r/30_refit-vpd.R fit brms
# models; r/33_simulate-null.R runs a 1,000-replicate Monte Carlo
# simulation (not brms, but still slow); r/04_calc-r2.R,
# r/05_summarize-pars.R, r/06_compare-gsw.R, and r/07_plot-curves.R each
# iterate over all ~2,100 individual weibull curve fits in
# objects/weibull/ (over a GB on disk) -- r/07_plot-curves.R additionally
# calls posterior_epred() per curve, making it the heaviest of the four --
# which is slow I/O/computation rather than a brms fit itself, but still
# too slow for `make fast`.

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
# NOTE: the dependency on $(STAMPS)/01_join-data is order-only (after the
# `|`), not a normal prerequisite. In the public repo, copy_files.R ships
# this stamp already future-dated so `make fast` can skip re-fitting; if
# 01_join-data were a normal prerequisite, Make would force this target to
# rebuild whenever 01_join-data itself needs rebuilding (which it always
# does in the copied repo, since joined-data.rds is intentionally
# regenerated fresh), regardless of this stamp's own timestamp. Order-only
# preserves correct build order for a from-scratch `make all` without
# defeating that skip mechanism. The same applies to the other brms/
# brms-adjacent targets below (03, 04, 05, 06, 10).
$(STAMPS)/02_fit-weibull: r/02_fit-weibull.R r/header.R r/functions.R \
  | $(STAMPS)/01_join-data $(STAMPS)
	@echo "==> r/02_fit-weibull.R (slow)"
	Rscript r/02_fit-weibull.R
	touch $@

# ---- brms: refit non-converged curves (slow) ----
$(STAMPS)/03_refit-weibull: r/03_refit-weibull.R r/header.R r/functions.R \
  | $(STAMPS)/02_fit-weibull $(STAMPS)
	@echo "==> r/03_refit-weibull.R (slow)"
	Rscript r/03_refit-weibull.R
	touch $@

$(STAMPS)/04_calc-r2: r/04_calc-r2.R r/header.R r/functions.R \
  | $(STAMPS)/03_refit-weibull $(STAMPS)
	@echo "==> r/04_calc-r2.R"
	Rscript r/04_calc-r2.R
	touch $@

$(STAMPS)/05_summarize-pars: r/05_summarize-pars.R r/header.R r/functions.R \
  | $(STAMPS)/03_refit-weibull $(STAMPS)
	@echo "==> r/05_summarize-pars.R"
	Rscript r/05_summarize-pars.R
	touch $@

$(STAMPS)/06_compare-gsw: r/06_compare-gsw.R r/header.R r/functions.R \
  | $(STAMPS)/03_refit-weibull $(STAMPS)
	@echo "==> r/06_compare-gsw.R"
	Rscript r/06_compare-gsw.R
	touch $@

# NOTE: $(STAMPS)/03_refit-weibull and $(STAMPS)/04_calc-r2 are order-only.
# A future-dated skip stamp (see copy_files.R) is always "newer" than any
# real timestamp, so if they were regular prerequisites this target would
# be considered out of date -- and thus rebuilt -- on every single `make
# fast` invocation, not just the first.
$(STAMPS)/07_plot-curves: r/07_plot-curves.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data | $(STAMPS)/03_refit-weibull $(STAMPS)/04_calc-r2 $(STAMPS)
	@echo "==> r/07_plot-curves.R"
	Rscript r/07_plot-curves.R
	touch $@

# NOTE: $(STAMPS)/05_summarize-pars is order-only; see the NOTE above
# r/07_plot-curves.R's rule.
$(STAMPS)/08_join-summary: r/08_join-summary.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data $(STAMPS)/01_join-data | $(STAMPS)/05_summarize-pars $(STAMPS)
	@echo "==> r/08_join-summary.R"
	Rscript r/08_join-summary.R
	touch $@

$(STAMPS)/09_make-tbl-vpd: r/09_make-tbl-vpd.R r/header.R r/functions.R \
  $(STAMPS)/00_load-data | $(STAMPS)
	@echo "==> r/09_make-tbl-vpd.R"
	Rscript r/09_make-tbl-vpd.R
	touch $@

# ---- brms: multiresponse phylogenetic model (slow) ----
# See the NOTE above r/02_fit-weibull.R's rule: these two prerequisites are
# order-only so the future-dated skip stamp shipped by copy_files.R works.
$(STAMPS)/10_fit-all: r/10_fit-all.R r/header.R r/functions.R \
  | $(STAMPS)/00_load-data $(STAMPS)/08_join-summary $(STAMPS)
	@echo "==> r/10_fit-all.R (slow)"
	Rscript r/10_fit-all.R
	touch $@

# NOTE: $(STAMPS)/10_fit-all is order-only; see the NOTE above
# r/07_plot-curves.R's rule.
$(STAMPS)/11_compare-models: r/11_compare-models.R r/header.R r/functions.R \
  | $(STAMPS)/10_fit-all $(STAMPS)
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

# NOTE: $(STAMPS)/10_fit-all is order-only; see the NOTE above
# r/07_plot-curves.R's rule.
$(STAMPS)/22_plot-collinear: r/22_plot-collinear.R r/header.R r/functions.R \
  | $(STAMPS)/10_fit-all $(STAMPS)
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

# ---- brms: refit with VPD covariate (slow) ----
# See the NOTE above r/02_fit-weibull.R's rule: these two prerequisites are
# order-only so a future-dated skip stamp (see copy_files.R) works, and
# `make fast` will use the pre-shipped objects/selected_model_vpd.rds
# instead of re-fitting.
$(STAMPS)/30_refit-vpd: r/30_refit-vpd.R r/header.R r/functions.R \
  | $(STAMPS)/12_select-model $(STAMPS)/27_summarize-vpd $(STAMPS)
	@echo "==> r/30_refit-vpd.R (slow)"
	Rscript r/30_refit-vpd.R
	touch $@

# NOTE: $(STAMPS)/30_refit-vpd is order-only; see the NOTE above
# r/07_plot-curves.R's rule.
$(STAMPS)/31_plot-tau-vpd-rate: r/31_plot-tau-vpd-rate.R r/header.R r/functions.R \
  $(STAMPS)/27_summarize-vpd | $(STAMPS)/30_refit-vpd $(STAMPS)
	@echo "==> r/31_plot-tau-vpd-rate.R"
	Rscript r/31_plot-tau-vpd-rate.R
	touch $@

# NOTE: $(STAMPS)/30_refit-vpd is order-only; see the NOTE above
# r/07_plot-curves.R's rule.
$(STAMPS)/32_plot-comparison: r/32_plot-comparison.R r/header.R r/functions.R \
  $(STAMPS)/12_select-model | $(STAMPS)/30_refit-vpd $(STAMPS)
	@echo "==> r/32_plot-comparison.R"
	Rscript r/32_plot-comparison.R
	touch $@

# ---- 1,000-replicate Monte Carlo null simulation (slow) ----
# See the NOTE above r/02_fit-weibull.R's rule: both prerequisites are
# order-only so a future-dated skip stamp (see copy_files.R) works, and
# `make fast` will use the pre-shipped objects/null-sim-fgmax-tau.rds
# instead of re-simulating.
$(STAMPS)/33_simulate-null: r/33_simulate-null.R r/header.R r/functions.R \
  | $(STAMPS)/01_join-data $(STAMPS)/05_summarize-pars $(STAMPS)
	@echo "==> r/33_simulate-null.R (slow)"
	Rscript r/33_simulate-null.R
	touch $@

# NOTE: $(STAMPS)/33_simulate-null is order-only; see the NOTE above
# r/07_plot-curves.R's rule.
$(STAMPS)/34_plot-null: r/34_plot-null.R r/header.R r/functions.R \
  $(STAMPS)/08_join-summary | $(STAMPS)/33_simulate-null $(STAMPS)
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

# All fast-script stamps, in pipeline order (used by both `fast` and `all`)
FAST_STAMPS := \
  $(STAMPS)/00_load-data \
  $(STAMPS)/01_join-data \
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
  $(STAMPS)/31_plot-tau-vpd-rate \
  $(STAMPS)/32_plot-comparison \
  $(STAMPS)/34_plot-null \
  $(STAMPS)/35_validate-nls-vs-bayes \
  $(STAMPS)/36_mediation-sensitivity \
  $(STAMPS)/37_plot-mediation-dag \
  $(STAMPS)/38_get-loo-diagnostics

# Slow scripts skipped by `make fast` (see the NOTE at the top of this
# file); not all of these are brms fits (r/33_simulate-null.R is a Monte
# Carlo simulation, and r/04_calc-r2.R/r/05_summarize-pars.R/
# r/06_compare-gsw.R/r/07_plot-curves.R are slow I/O and/or posterior
# computation over objects/weibull/ rather than a model fit), but all are
# slow enough to pre-ship their output.
SLOW_STAMPS := \
  $(STAMPS)/02_fit-weibull \
  $(STAMPS)/03_refit-weibull \
  $(STAMPS)/04_calc-r2 \
  $(STAMPS)/05_summarize-pars \
  $(STAMPS)/06_compare-gsw \
  $(STAMPS)/07_plot-curves \
  $(STAMPS)/10_fit-all \
  $(STAMPS)/30_refit-vpd \
  $(STAMPS)/33_simulate-null

# --------------------------------------------------------------------------
# Top-level targets
# --------------------------------------------------------------------------

# Render ms/ms.pdf using whatever outputs already exist -- no R scripts run
pdf:
	$(RENDER)

# Run all fast scripts (requires the slow scripts' outputs to already
# exist from a prior `make all`; see the NOTE at the top of this file),
# then render
fast: $(FAST_STAMPS)
	$(RENDER)

# Run every script (r/00_ through r/38_) including the slow ones, then render
all: $(SLOW_STAMPS) $(FAST_STAMPS)
	$(RENDER)

# --------------------------------------------------------------------------
# Diff: generate a tracked-changes PDF comparing submitted vs revised ms.
# Requires latexdiff (bundled with most TeX distributions).
# Run `make pdf` first so ms/ms.tex is up to date, then run `make diff`.
# --------------------------------------------------------------------------

.PHONY: diff
diff: ms/ms-v1-submitted.tex ms/ms.tex
	latexdiff --encoding=utf8 \
	  --add-to-config "PICTUREENV=longtable[\w\d*@]*" \
	  ms/ms-v1-submitted.tex ms/ms.tex > ms/ms-diff.tex
	cd ms && latexmk -pdflua -interaction=nonstopmode ms-diff.tex
	cd ms && latexmk -c ms-diff.tex

# --------------------------------------------------------------------------
clean:
	rm -rf $(STAMPS)
	rm -f ms/ms.pdf
	rm -f ms/ms-diff.tex ms/ms-diff.pdf
