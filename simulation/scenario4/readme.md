
Run

```r
Rscript .\scenario4_boundary_oracle.R `
  --n=200,500,1000,2000,5000,10000,20000 `
  --B=500 `
  --cores=4 `
  --c-h=0.8 `
  --tail-p=0 `
  --grid-max=262144 `
  --profile-radius-se=8 `
  --profile-grid-points=65 `
  --qq-n=20000 `
  --out-dir='scenario4_boundary_oracle_resultsB500' `
  --cache-dir='oracle_nuisance_fixed_resultsB500'
```
