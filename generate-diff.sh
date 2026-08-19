#!/usr/bin/env bash
# Generate ms/ms-diff.pdf: a tracked-changes PDF comparing the originally
# submitted manuscript (ms/ms-v1-submitted.tex) against the current
# revision (rendered fresh from ms/ms.qmd).
#
# This is a thin wrapper around the `pdf` and `diff` Makefile targets
# (kept here for convenience running outside of `make`, e.g. as a
# standalone step), so the actual latexdiff/latexmk logic has a single
# source of truth in the Makefile.
#
# Usage: ./generate-diff.sh
# Requires: quarto, GNU make, latexdiff, and a lualatex-capable latexmk
# (this manuscript needs lualatex specifically, for mainfont: Times New
# Roman via fontspec and non-ASCII characters like the 'okina in
# "Hawai'i").

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

if [ ! -f ms/ms-v1-submitted.tex ]; then
  echo "ERROR: ms/ms-v1-submitted.tex not found; nothing to diff against." >&2
  exit 1
fi

echo "==> Rendering current manuscript (make pdf)"
make pdf

echo "==> Building tracked-changes diff (make diff)"
make diff

if [ -f ms/ms-diff.pdf ]; then
  echo "==> Done: ms/ms-diff.pdf"
else
  echo "ERROR: ms/ms-diff.pdf was not created; see ms/ms-diff.log" >&2
  exit 1
fi
