# Field Evaluation of Integrated Soybean Cyst Nematode Management Using Spatially Informed Mixed Models 

## Citation
The following source code accompanies our publication, to be submitted to *Phytopathology* Journal:

```
@article{garnica202x,
  author = {Vinicius C. Garnica and Horacio D. Lopez-Nicora},
  title = {Field evaluation of integrated soybean cyst nematode management using spatially informed mixed models},
  year = {202x},
  doi = {xxx},
  journal = {xxx}
}
```

## Introduction

This repository contains data and code for analyzing the effectiveness of integrated management 
strategies for soybean cyst nematode (SCN) using spatially informed mixed models. The study evaluates 
SCN population dynamics (reproduction factor, *rf*) and yield responses (*yld*) across field trials 
conducted in 2022 and 2023 in Ohio, incorporating spatial autocorrelation structures to 
account for within-field heterogeneity and improve estimation precision.

The analysis employs advanced mixed-model frameworks to:

* Model spatial dependencies in field data using various covariance structures
* Evaluate the performance of SCN management strategies including resistant cultivars and seed treatments
* Account for plot-level spatial trends and environmental heterogeneity
* Quantify cultivar-specific damage coefficients relating initial SCN population to yield loss
* Provide robust estimates of treatment effects adjusted for spatial autocorrelation

This repository outlines the complete analysis pipeline, detailing the necessary Quarto markdown documents, R scripts, and data files for reproducing the study. The folder structure is organized as follows:

```
SCN_spatial_inference/
├── spatial_modeling_rf.qmd              # Analysis of SCN reproduction factor
├── spatial_modeling_rf.html             # Rendered HTML output for rf analysis
├── spatial_modeling_yld.qmd             # Analysis of soybean yield responses
├── spatial_modeling_yld.html            # Rendered HTML output for yield analysis
├── pioneers.csv                         # Dataset
├── functions.R                          # Custom functions for spatial modeling
├── results/
│     ├── cult_int/                      # Cultivar and cultivar-by-seed treatment interaction estimates
│     │     ├─── cult_int_2022.csv
│     │     ├─── cult_int_2022_rf.csv
│     │     ├─── cult_int_2023.csv
│     │     ├─── cult_int_2023_rf.csv
│     │     ├─── cult_int_multi.csv      # Multi-year analysis
│     │     └─── cult_int_multi_rf.csv
│     ├── damage/                        # Damage function coefficients
│     │     ├─── coef_damage_2022.csv
│     │     ├─── coef_damage_2023.csv
│     │     └─── coef_damage_multi.csv
│     ├── model_selection/               # Spatial model comparison
│     │     ├─── comparison_2022.csv
│     │     ├─── comparison_2022_rf.csv
│     │     ├─── comparison_2023.csv
│     │     ├─── comparison_2023_rf.csv
│     │     ├─── comparison_multi.csv
│     │     └─── comparison_multi_rf.csv
│     ├── seed_trt/                      # Seed treatment effect estimates
│     │     ├─── seed_trt_2022.csv
│     │     ├─── seed_trt_2022_rf.csv
│     │     ├─── seed_trt_2023.csv
│     │     ├─── seed_trt_2023_rf.csv
│     │     ├─── seed_trt_multi.csv
│     │     └─── seed_trt_multi_rf.csv
│     ├── variances/                     # Variance component estimates
│     │     ├─── var_table_2022.csv
│     │     ├─── var_table_2022_rf.csv
│     │     ├─── var_table_2023.csv
│     │     ├─── var_table_2023_rf.csv
│     │     ├─── var_table_multi.csv
│     │     └─── var_table_multi_rf.csv
│     ├── drawing.svg                    # Conceptual diagram
│     ├── fig_1.tiff
│     ├── fig_2.tiff
│     ├── fig_3.tiff
│     ├── fig_4.tiff
│     └── fig_S1.tiff
└── functions.R

```
## Analysis Pipeline

### Yield Response Analysis (`spatial_modeling_yld.qmd`)

This Quarto document analyzes soybean yield data from 2022, 2023, and combined multi-year trials. The 
workflow imports and prepares yield data, fits mixed models with different spatial structures, and 
exports CSV files with model comparison metrics (AIC/BIC), treatment effect estimates (saved in `cult_int/` and `seed_trt/`), and 
variance components (saved in `variances/`).

### Reproduction Factor Analysis (`spatial_modeling_rf.qmd`)

This Quarto document analyzes the SCN reproduction factor (rf = log(Pf/Pi)) using the same modeling framework as the yield analysis.

### Custom Functions (`functions.R`)

This code defines a suite of R functions to evaluate ASReml models and extract biologically interpretable contrasts. The `icREML` function 
computes AIC and BIC from a list of fitted models using restricted maximum likelihood. The `st_contrasts`, `cult_contrasts_ss`, and `cult_contrasts_ms` 
functions estimate and optionally back-transform contrasts for seed treatment and cultivar effects, while `spatial_matrix` constructs 
spatial basis matrices for tensor-product spline modeling across field trial sites.

## References

* Gilmour, A. R. and Cullis, B. R. (1997). Accounting for natural and extraneous variation in the analysis of field experiments. Journal of Agricultural, Biological, and Environmental Statistics, 2, 269–293.

* Velazco, J. G., Rodríguez-Álvarez, M. X., Boer, M. P., Jordan, D. R., Eilers, P. H. C., Malosetti, M., and van Eeuwijk, F. A. (2017). Modelling spatial trends in sorghum breeding field trials using a two-dimensional P-spline mixed model. Theoretical and Applied Genetics, 130, 1375–1392.

* Welham, S. J. (2022). TPSbits: Creates structures to enable fitting and examination of 2D tensor-product splines using ASReml-R (R package version 1.0.2) [Computer software]. Retrieved August 4, 2025, from https://mmade.org/tpsbits
