# Galaxies Final Analysis

This folder contains the final empirical analysis script for the galaxies data:

```r
galaxies_final_analysis.R
```

The script analyzes `MASS::galaxies` using PMHD-MCP order selection, PMHD-MCP-R fixed-order refit, and a Gaussian-mixture comparison.

## Requirements

R packages:

```r
install.packages(c("MASS", "nloptr", "mclust"))
```

## Default Run

From this folder:

```bash
Rscript galaxies_final_analysis.R
```

By default, the script uses the locked paper results and writes them to:

```text
K6final/
```

## Fast Run

```bash
GALAXIES_FAST_MODE=1 Rscript galaxies_final_analysis.R
```

On Windows PowerShell:

```powershell
$env:GALAXIES_FAST_MODE="1"
Rscript .\galaxies_final_analysis.R
```

This writes the same locked final results to:

```text
K6fast/
```

## Recompute Instead of Copying Locked Results

In R:

```r
source("galaxies_final_analysis.R")

run_galaxies_final_analysis(
  output_dir = "K6final",
  use_locked_results = FALSE,
  force_path_refit = TRUE
)
```

For a faster fixed-\(K=6\), fixed-\(\lambda\) warm-start recomputation:

```r
run_galaxies_final_analysis(
  output_dir = "K6fast",
  use_locked_results = FALSE,
  fast_paper_mode = TRUE,
  force_path_refit = TRUE
)
```

## Main Outputs

- `galaxies_pmhd_refit_summary.csv`
- `galaxies_pmhd_refit_table.csv`
- `galaxies_pmhd_refit_fit.pdf`
- `galaxies_selector_paths.pdf`
- `galaxies_pmhd_refit_results.rds`

The paper reports the locked final results:

```text
PMHD-MCP raw   K=6  lambda=0.01564243  H=0.10409255
PMHD-MCP-R     K=6                    H=0.02322901
Gaussian mix   K=4                    H=0.24255342
```
