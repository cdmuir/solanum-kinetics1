# --- Download and preprocess raw data from the companion solanum-aa repository ---
#
# Each file below is pinned to a specific solanum-aa commit (not `main`) so
# that re-running this script fetches the same upstream data every time,
# rather than whatever happens to be on solanum-aa's main branch that day.
#
# IMPORTANT: the pins are NOT all the same commit. Before this fix, this
# script always read from `main`, and a `git log` on each of these files in
# *this* repo showed they were actually last re-fetched/committed at
# different times (phylogeny.rds: 2025-10-10; accession_info.rds/
# plant_info.rds: 2025-10-21; df_germ_summary.rds/df_growth_summary.rds:
# 2026-02-19; rh_curves.rds/stomata.rds: 2026-03-18) -- i.e. these files
# were never a single upstream snapshot to begin with. Each pin below is
# the solanum-aa commit that was HEAD for that specific file's path as of
# its corresponding date here, found via
# `commits?path=<file>&until=<date>`. Only phylogeny.rds's pin was
# verified as an exact byte-for-byte match (via sha256, against the
# existing Git LFS pointer in this repo); the rest could not be verified
# that way because this script re-serializes every file through
# write_rds() (and transforms a few via mutate()/select()), so a raw-byte
# hash comparison isn't meaningful for them -- their *logical* content
# should still match, since the pin was chosen from this repo's own commit
# history for that exact file, not guessed.
#
# To intentionally pick up newer upstream data for one or more files,
# update the relevant REF constant below (check
# https://github.com/cdmuir/solanum-aa/commits/main for candidates) and
# rerun `make all` from scratch.
#
# NOTE: solanum-aa's `v0.1` tag is NOT a safe alternative to any of these
# pins -- it is far behind main and predates several data bug fixes (e.g.
# "fixed bug in stomata.rds", "corrected stomatal density for edge
# effect", "filtered data with bad Tleaf").

source("r/header.R")

solanum_aa_url <- function(ref, path) {
  glue("https://github.com/cdmuir/solanum-aa/raw/{ref}/{path}")
}

# accession-info.rds / plant-info.rds: as of 2025-10-21 ("updating data")
ACCESSION_PLANT_INFO_REF <- "a1b5a64d2ebe1f734cb66b38ad68d70c5025162c"

read_rds(solanum_aa_url(ACCESSION_PLANT_INFO_REF, "data/accession-info.rds")) |>
  write_rds("data/accession_info.rds")

read_rds(solanum_aa_url(ACCESSION_PLANT_INFO_REF, "data/plant-info.rds")) |>
  write_rds("data/plant_info.rds")

# phylogeny.rds: as of 2025-10-10 ("working on analyzing tau") -- verified
# exact sha256 match against this repo's committed Git LFS pointer
PHYLOGENY_REF <- "f8988102604adb49c4ceb56fd8319bfbf1c4edf6"

read_rds(solanum_aa_url(PHYLOGENY_REF, "data/phylogeny.rds")) |>
  write_rds("data/phylogeny.rds")

# trimmed_rh_curves.rds / stomata.rds: as of 2026-03-18 ("updating fit-gcl")
RH_STOMATA_REF <- "3a37cd54ac3258871c7f225479a73b4fce68fc68"

read_rds(solanum_aa_url(RH_STOMATA_REF, "data/trimmed_rh_curves.rds")) |>
  mutate(ci = as.numeric(as.factor(curve))) |>
  mutate(t_sec = elapsed - min(elapsed), .by = ci) |>
  write_rds("data/rh_curves.rds")

read_rds(solanum_aa_url(RH_STOMATA_REF, "data/stomata.rds")) |>
  select(-contains("pavement"), -contains("index")) |>
  write_rds("data/stomata.rds")

# df_germ_summary.rds / df_growth_summary.rds: as of 2026-02-19 ("updating methods")
GERM_GROWTH_SUMMARY_REF <- "4ede1ac282da270d92f2c7835bb482fa0d5a43aa"

read_rds(solanum_aa_url(GERM_GROWTH_SUMMARY_REF, "data/df_germ_summary.rds")) |>
  write_rds("data/df_germ_summary.rds")

read_rds(solanum_aa_url(GERM_GROWTH_SUMMARY_REF, "data/df_growth_summary.rds")) |>
  write_rds("data/df_growth_summary.rds")
