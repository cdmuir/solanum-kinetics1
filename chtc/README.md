# Running `r/10_fit-all.R`'s 64 models on CHTC

This folder lets you run each of the 64 brms models fit locally by
`r/10_fit-all.R` as an independent HTCondor job on UW–Madison's CHTC HTC
system, each requesting 3 CPUs so its 3 chains sample in parallel
(`cores = 3, chains = 3`).

`r/10_fit-all.R` itself is unchanged and still works for local runs;
everything here is a self-contained parallel workflow.

## Files

| File | Purpose |
|---|---|
| `image.def` | Apptainer definition: R + brms + cmdstanr + a compiled CmdStan, built once |
| `fit_model.R` | Fits ONE model, given an index; the job executable |
| `test_one.sub` | Submit file for a single test job (run this first!) |
| `fit_all.sub` | Submit file for all 64 jobs |
| `reassemble_results.R` | Run locally after downloading results; rebuilds `objects/fits/` + `objects/df_forms.rds` |

## Prerequisites

- A CHTC account and access to a CHTC access point (submit node), connected
  via UW network/VPN with 2FA. See
  [Log into CHTC](https://chtc.cs.wisc.edu/uw-research-computing/connecting).
- A `/staging` directory on CHTC (request one from `chtc@cs.wisc.edu` if you
  don't have one — needed for the container image).

## Step 1: Build and test the container

Containers must be built on a CHTC build node via an interactive job.

1. Copy `chtc/image.def` to your CHTC access point (e.g. via `scp`).
2. Start an interactive build job:
   ```
   # build.sub
   log = build.log
   transfer_input_files = image.def
   +IsBuildJob = true
   request_cpus = 4
   request_memory = 8GB
   request_disk = 10GB
   queue
   ```
   ```
   condor_submit -i build.sub
   ```
3. Inside the interactive job:
   ```
   apptainer build solanum-fit.sif image.def
   ```
   This installs R, brms, cmdstanr, and compiles CmdStan — expect this to
   take a while (mostly CmdStan compilation).
4. Test it:
   ```
   apptainer shell -e solanum-fit.sif
   R -e 'library(brms); library(cmdstanr); cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN")); cmdstanr::cmdstan_version()'
   exit
   ```
5. Move the image to staging and exit the interactive job:
   ```
   mv solanum-fit.sif /staging/$USER/
   exit
   ```

## Step 2: Upload scripts and data

From your local machine, copy the following to a working directory on the
CHTC access point (e.g. `~/solanum-kinetics1-chtc/`):

- `chtc/fit_model.R`
- `chtc/fit_all.sub`
- `chtc/test_one.sub`
- `data/joined-summary.rds`
- `data/phylogeny.rds`

Create a `logs/` subdirectory there too (`mkdir logs`).

Edit `test_one.sub` and `fit_all.sub`: replace `<user>` in the
`container_image` path with your CHTC username (or wherever you placed
`solanum-fit.sif` under `/staging`).

## Step 3: Test with ONE job first

**Do not skip this.** Resource requests in `fit_all.sub` (16GB memory, 12GB
disk) are unverified placeholders.

```
condor_submit test_one.sub
condor_q                      # watch it run
```

When it finishes, check:

```
cat logs/test_0.log           # look for "Memory (MB)" and "Disk (KB)" usage
cat logs/test_0.err           # check for errors/warnings
ls -la fit_01.rds             # confirm the output came back
```

If the job held or failed, run `condor_q -hold` to see why, and check
`logs/test_0.err` for R/Stan errors.

Update `request_memory` / `request_disk` in `fit_all.sub` based on the
actual usage reported in `logs/test_0.log` (with some margin), before
proceeding.

## Step 4: Submit all 64 jobs

```
condor_submit fit_all.sub
condor_q
```

Each job writes `fit_<NN>.rds` (`fit_01.rds`–`fit_64.rds`) to the working
directory, which HTCondor transfers back automatically on completion.

Monitor with `condor_q`, and if any jobs go on hold, inspect them with
`condor_q -hold` and the corresponding `logs/fit_<N>.log`/`.err` files.

## Step 5: Retrieve results and reassemble locally

Once all 64 jobs complete, copy the `fit_*.rds` files back to your local
machine, e.g.:

```
# from your local machine
mkdir -p chtc/downloaded_fits
scp '<user>@ap2001.chtc.wisc.edu:~/solanum-kinetics1-chtc/fit_*.rds' chtc/downloaded_fits/
```

Then, from the project root locally:

```r
# in R, or:
# Rscript chtc/reassemble_results.R chtc/downloaded_fits
```

This validates all 64 files load as `brmsfit` objects, copies them into
`objects/fits/`, and rewrites `objects/df_forms.rds` using the same
`(lambda_idx, tau_idx)` mapping as `chtc/fit_model.R`, so
`r/11_compare-models.R` and downstream scripts work unchanged.

**Note:** this overwrites `objects/fits/` and `objects/df_forms.rds`, which
currently hold a prior local run. Back them up first if you want to keep
both.

## Index mapping reference

Both `fit_model.R` and `reassemble_results.R` use the same explicit mapping
from a 0-based HTCondor `$(Process)` (0–63) to a (lambda formula, tau
formula) pair, rather than relying on `tidyr::crossing()`'s ordering:

```r
i <- Process + 1          # 1..64
lambda_idx <- (i - 1) %/% 8   # 0..7 -> bf_lambda0..bf_lambda7
tau_idx    <- (i - 1) %% 8    # 0..7 -> bf_tau0..bf_tau7
```

`fit_01.rds` = `bf_lambda0` + `bf_tau0`, `fit_02.rds` = `bf_lambda0` +
`bf_tau1`, ..., `fit_09.rds` = `bf_lambda1` + `bf_tau0`, etc.
