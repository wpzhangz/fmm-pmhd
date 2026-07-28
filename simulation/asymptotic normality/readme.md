# Ensure the core library is available (it requires simulation_stability_elbow_parallel_main.R)
```r
ls simulation_stability_elbow_parallel_main.R
```
# Syntax and loading check
```r
Rscript -e 'source("scenario3_asymptotic_normality_final.R"); cat("loaded OK\n")'
```
# Smoke test (single core, 20 iterations, n=500)
```r
Rscript scenario3_asymptotic_normality_final.R --quick
```
# Formal run (should be consistent with previous s1d results)
```r
Rscript scenario1c_asymptotic_normality_final.R \
  --scenario=s1d --n=500,1000,2000 --B=200 --cores=12 --bw-constant=2.0
python3 make_normality_table_s1d.py normality_s1d --bwc=2.0
```
