#!/usr/bin/env Rscript

# Script to assemble the Dryad data package into ./dryad
#
# Per the manuscript's Data Availability statement (ms/ms.qmd), the Dryad
# archive consists of the two curve- and accession-level estimate tables
# referenced in @tbl-estimates-curve and @tbl-estimates-accession, each with
# its accompanying data dictionary, plus a README describing the dataset.
# Raw phenotype data are hosted separately in cdmuir/solanum-aa and are not
# duplicated here; analysis code is archived separately on Zenodo/GitHub.

dest_dir <- "dryad"

# Create destination directory if it doesn't exist
if (!dir.exists(dest_dir)) {
  dir.create(dest_dir, recursive = TRUE)
  cat("Created directory:", dest_dir, "\n")
}

# Clear out any existing contents of dest_dir before copying fresh files, so
# stale files from a previous run don't linger.
existing_items <- list.files(dest_dir, all.files = TRUE, no.. = TRUE, full.names = FALSE)
if (length(existing_items) > 0) {
  unlink(file.path(dest_dir, existing_items), recursive = TRUE, force = TRUE)
  cat("Cleared", length(existing_items), "existing item(s) from", dest_dir, "\n\n")
}

# Data files and their dictionaries, as described in ms/ms.qmd
# (@tbl-estimates-curve, @tbl-estimates-accession)
files_to_copy <- c(
  "tables/tbl-estimates-curve.csv",
  "tables/tbl-estimates-curve-dictionary.csv",
  "tables/tbl-estimates-accession.csv",
  "tables/tbl-estimates-accession-dictionary.csv"
)

cat("Copying data files to", dest_dir, "...\n\n")

success_count <- 0
for (file in files_to_copy) {
  src_path <- file.path(getwd(), file)
  dest_path <- file.path(dest_dir, basename(file))

  if (!file.exists(src_path)) {
    cat("WARNING: Source not found:", src_path, "\n")
    next
  }

  file.copy(src_path, dest_path, overwrite = TRUE)
  cat("Copied file:", file, "\n")
  success_count <- success_count + 1
}

cat("\n\u2713 Successfully copied", success_count, "of", length(files_to_copy), "data files to", dest_dir, "\n\n")

# Write the Dryad README describing the dataset
readme_text <- paste0(
"# Data for: Guard cell size and initial conductance influence stomatal closure kinetics\n\n",
"Christopher D. Muir, Wei Shen Lim\n\n",
"This dataset contains the derived data tables underlying the results reported in the manuscript. ",
"Raw phenotype data (stomatal anatomy, growth conditions, and germination records) are hosted separately in the ",
"[cdmuir/solanum-aa](https://github.com/cdmuir/solanum-aa) repository and remain subject to that repository's own terms. ",
"Analysis code is available at [github.com/cdmuir/solanum-kinetics](https://github.com/cdmuir/solanum-kinetics) and archived on Zenodo upon publication.\n\n",
"## Description of the data and file structure\n\n",
"| File | Description |\n",
"|------|-------------|\n",
"| `tbl-estimates-curve.csv` | Estimates of stomatal anatomy and kinetic parameters (guard cell length, initial and maximum stomatal conductance, lag time, and time constant) associated with each individual humidity response curve. These are the data used to fit the multiresponse models reported in the manuscript. |\n",
"| `tbl-estimates-curve-dictionary.csv` | Data dictionary describing each column in `tbl-estimates-curve.csv`. |\n",
"| `tbl-estimates-accession.csv` | Population (accession)-level estimates of stomatal anatomy and kinetic parameters, based on multiresponse model predictions. |\n",
"| `tbl-estimates-accession-dictionary.csv` | Data dictionary describing each column in `tbl-estimates-accession.csv`. |\n\n",
"Each data dictionary lists, for every column in the corresponding CSV, the variable name, data type, acceptable values (for categorical variables), and a description.\n\n",
"## Sharing/access information\n\n",
"- License: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)\n",
"- Related publication: Muir CD, WS Lim. Guard cell size and initial conductance influence stomatal closure kinetics.\n",
"- Code to reproduce all analyses: [github.com/cdmuir/solanum-kinetics](https://github.com/cdmuir/solanum-kinetics)\n"
)

writeLines(readme_text, file.path(dest_dir, "README.md"))
cat("Wrote README.md to", dest_dir, "\n")
