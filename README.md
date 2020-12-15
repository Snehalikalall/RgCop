# RgCop
Regularized Copula based Feature Selection- Application on Single cell RNA Sequence Data


Load the dataset from **Data** folder, and put the codes and Datasets in same path.
Run the DataProcessing.R code. 
Run the FeatureSelection.R in **R** folder. 

The default core, p=20, and Tuning Coefficient-value is (coeff=0.009) for RgCop. 

The MarkerSelection_RgCop.ipynb is for generating all the clustering validation figures.

## Pre-requisites

> R version  4.0.2

> Python 3.7

> R packages: SingleCellExperiment, edgeR, scDatasets, biomaRt, Linnorm, Rfast, foreach, doParallel, MASS

> Python Packages: scanpy, leidenalg



