## High-Dimensional Locally Stationary Models for Identifying Neuroimaging Biomarkers in Autism Spectrum Disorder

### Introduction

This page includes data and functions associated with the paper "High-dimensional locally stationary models for identifying neuroimaging biomarkers in Autism Spectrum Disorder". 

#### Data

`Caltech_0051475_rois_aal.1D` is a rs-fMRI dataset for a control subject from the Autism Brain Image Exchange (ABIDE) repository. For more information about this repository, please visit [ABIDE Preprocessed](preprocessed-connectomes-project.org/abide/).

#### Functions

`h2part` can be used to conduct hypothesis test comparing distributions (Poisson, Negative Binomial, Zero-Inflated Poisson, and Zero-Inflated Negative Binomial) of edge counts for ASD group and control group according to paper Togo et al.

`hboot` can be used to conduct hypothesis test comparing distributions (Poisson, Negative Binomial, Zero-Inflated Poisson, and Zero-Inflated Negative Binomial) of edge counts for ASD group and control group using bootstrap confidence intervals. 

`hboot2` can be used to conduct hypothesis test comparing distributions (Poisson and Negative Binomial) of edge counts for ASD group and control group using bootstrap confidence intervals. 

`mconn` can be used to calculate connectivity matrix across all frequencies. 

`mconn2` can be used to calculate connectivity matrix across certain frequencies. 

`mnet` can be used to calculate number of edges for all brain networks. 

`medge` can be used to calculate number of edges for certain brain regions. 





