#!/usr/bin/env Rscript

# Script to copy all necessary files to ../solanum-kinetics directory
# Based on dependencies listed in Makefile

# Destination directory
dest_dir <- "../solanum-kinetics"

# Create destination directory if it doesn't exist
if (!dir.exists(dest_dir)) {
  dir.create(dest_dir, recursive = TRUE)
  cat("Created directory:", dest_dir, "\n")
}

# Define all files needed based on Makefile dependencies
files_to_copy <- c(
  # R scripts
  "r/00_load-data.R",
  "r/01_join-data.R",
  "r/02_fit-weibull.R",
  "r/03_refit-weibull.R",
  "r/04_calc-r2.R",
  "r/05_summarize-pars.R",
  "r/06_compare-gsw.R",
  "r/07_plot-curves.R",
  "r/08_join-summary.R",
  "r/09_make-tbl-vpd.R",
  "r/10_fit-all.R",
  "r/11_compare-models.R",
  "r/12_get-partial-cor.R",
  "r/13_make-tbl-estimates-curve.R",
  "r/14_make-tbl-estimates-accession.R",
  "r/15_plot-accession-anatomy.R",
  "r/16_make-tbl-fit-summary.R",
  "r/17_plot-accession-kinetics.R",
  "r/18_plot-variance.R",
  "r/19_plot-gcl-tau.R",
  "r/20_plot-colinear.R",
  "r/21_plot-fgmax-kinetics.R",
  "r/22_plot-mediation.R",
  "r/23_plot-gcl.R",
  "r/24_make-illustrations.R",
  "r/25_plot-conceptual.R",
  "r/header.R",
  "r/functions.R",
  
  # Manuscript
  "ms/ms.qmd",
  
  # Data and pre-computed outputs
  "data/",
  "objects/",
  
  # Git configuration files
  ".gitignore",
  ".gitattributes",
  
  # Other files
  "Makefile",
  "README.md"
)

# Function to copy files and directories
copy_file_or_dir <- function(src, dest_base) {
  src_path <- file.path(getwd(), src)
  dest_path <- file.path(dest_base, src)

  if (!file.exists(src_path)) {
    cat("WARNING: Source not found:", src_path, "\n")
    return(FALSE)
  }

  if (dir.exists(src_path)) {
    if (!dir.exists(dest_path)) {
      dir.create(dest_path, recursive = TRUE)
    }
    # Use rsync for directories; exclude objects/weibull/
    exclude <- if (src == "objects/") "--exclude='weibull'" else ""
    cmd <- paste("rsync -a", exclude, shQuote(src_path), shQuote(dest_path))
    system(cmd)
    cat("Copied directory:", src, "\n")
  } else {
    # Copy file, creating parent directory if needed
    dest_parent <- dirname(dest_path)
    if (!dir.exists(dest_parent)) {
      dir.create(dest_parent, recursive = TRUE)
    }
    file.copy(src_path, dest_path, overwrite = TRUE)
    cat("Copied file:", src, "\n")
  }

  return(TRUE)
}

# Copy all files
cat("Copying files to", dest_dir, "...\n\n")

success_count <- 0
for (file in files_to_copy) {
  if (copy_file_or_dir(file, dest_dir)) {
    success_count <- success_count + 1
  }
}

cat("\n✓ Successfully copied", success_count, "items to", dest_dir, "\n")
