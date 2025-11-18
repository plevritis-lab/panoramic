# PANORAMIC
**P**ooled **AN**alysis **O**f Va**R**iance-**A**ware **M**odeling and **I**nference of **C**olocalization (**PANORAMIC**)

This package is currently under development for submission to Bioconductor. 

![alt text](https://github.com/plevritis-lab/panoramic/blob/main/images/panoramic_graphic.png)

## Installation
```r
devtools::install_github("plevritis-lab/panoramic")
library(panoramic)
```

## Tutorial
Start by preparing the inputs to panoramic: (1) a list of SpatialExperiment objects containing your spatial data and (2) a design matrix mapping samples to group level labels for comparison. If no comparison is being performed, then assign all samples to a single group. 
```r
```

Next, create the panoramic SummarizedExperiment object compiled by calcualting relevant spatial statistics across all pairs of cell types
```r
se <- panoramic(
  spe_list,
  design,
  cell_type="cell_type",
  radii_um=c(75),
  stat="Lcross",
  nsim=100,
  BPPARAM=BiocParallel::MulticoreParam(parallel::detectCores())
)  
```

Synthesize spatial informations across samples in a group, then test for differences
```r
se_meta <- panoramic_meta(se, tau2="SJ")
se_diff <- panoramic_test_groups(se_meta)
```
