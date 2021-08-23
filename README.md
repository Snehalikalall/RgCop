# RgCop
Regularized Copula based Feature Selection- Application on Single cell RNA Sequence Data



## Pre-requisites

> R version  4.0.2

> Python 3.7

> Python packages: scanpy

> R packages: copula, foreach, doParallel

## Install
library("devtools")

install_github("Snehalikalall/RgCop")

Check the installation:

library(RgCop)


## Load required packages

R packages

     library(SingleCellExperiment)
     library(foreach)
     library(doParallel)
     library(Linnorm)
     library(copula)

Python Packages: 
 
    pip install scanpy



## Preprocessing of raw data

Preprocess the raw data using DataProcessing.R function

    Biase_data<- readRDS("muraro.rds")
    data <- assay(Biase_data) 
    annotation <- Biase_data[[1]] #already factor type class
    colnames(data) <- annotation
    data_process = normalized_data(data)
    write.table(t(data_process),file="RgCop/Data/data_process.csv",sep=",",row.names = FALSE,col.names = FALSE)


## Select genes feom sampled data 

Use RgCopfeature.R to select informative feature subset using Regularized copula based feature selection

    # load the "pollen_process.csv" and "pollenc.csv" in R file.  
    data=as.matrix(read.csv("Data/pollen_process.csv",header=FALSE))    #should be cells in row, genes in coloumn.
    classd<-as.matrix(read.csv("Data/pollenc.csv",header=FALSE))
    p=40
    nf=500
    #nf: Number of feature to be selected, default is 500; P: Number of cores, default is 40

    # Regularized copula based feature selection, the function returns two elements: i) Data with selected features, and ii) The selected feature subset
    Result=RgCopfeature(data,classd, p,nf)
    # Data with selected features
    Result$Feadata
    # The selected feature subset
    Result$Features
