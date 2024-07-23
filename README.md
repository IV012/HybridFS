## R package `HybridFS` 

### Description

This R package implements: 
  1. HFS, an efficient nonlinear feature screening method using multiple utilities.
  2. HiFIT, a high-dimensional feature importance testing via machine learning models.

### Minimal Requirement

- System requirement: Rtools (Windows), None (Linux, Mac)
  
- Dependent on R (>= 4.0.0), methods, e1071, glmnet, randomForest, xgboost, foreach, parallel, doParallel, dHSIC, isotree

- Multi-core CPU for parallel computing

### Installation

```{r}
devtools::install_github("IV012/HybridFS")
```

### Related Pages

- [Deep Treatment Learning](https://github.com/SkadiEye/deepTL/tree/master) by Xinlei Mi

### References

Mi, X., Zou, B., Zou, F. et al. Permutation-based identification of important biomarkers for complex diseases via machine learning models. Nat Commun 12, 3008 (2021). https://doi.org/10.1038/s41467-021-22756-2

