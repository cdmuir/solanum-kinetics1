# Data for: Guard cell size and initial conductance influence stomatal closure kinetics

Christopher D. Muir, Wei Shen Lim

This dataset contains the derived data tables underlying the results reported in the manuscript. Raw phenotype data (stomatal anatomy, growth conditions, and germination records) are hosted separately in the [cdmuir/solanum-aa](https://github.com/cdmuir/solanum-aa) repository and remain subject to that repository's own terms. Analysis code is available at [github.com/cdmuir/solanum-kinetics](https://github.com/cdmuir/solanum-kinetics) and archived at [10.5281/zenodo.22802259](https://doi.org/10.5281/zenodo.22802259).

## Description of the data and file structure

| File | Description |
|------|-------------|
| `tbl-estimates-curve.csv` | Estimates of stomatal anatomy and kinetic parameters (guard cell length, initial and maximum stomatal conductance, lag time, and time constant) associated with each individual humidity response curve. These are the data used to fit the multiresponse models reported in the manuscript. |
| `tbl-estimates-curve-dictionary.csv` | Data dictionary describing each column in `tbl-estimates-curve.csv`. |
| `tbl-estimates-accession.csv` | Population (accession)-level estimates of stomatal anatomy and kinetic parameters, based on multiresponse model predictions. |
| `tbl-estimates-accession-dictionary.csv` | Data dictionary describing each column in `tbl-estimates-accession.csv`. |

Each data dictionary lists, for every column in the corresponding CSV, the variable name, data type, acceptable values (for categorical variables), and a description.

## Sharing/access information

- License: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
- Related publication: Muir CD, WS Lim. Guard cell size and initial conductance influence stomatal closure kinetics.
- Code to reproduce all analyses: [github.com/cdmuir/solanum-kinetics](https://github.com/cdmuir/solanum-kinetics)

