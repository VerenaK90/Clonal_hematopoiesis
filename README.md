# Clonal_hemopoiesis
This repository contains code to reproduce analysis in Körber et al., Quantifying drift and selection from a single bulk tissue sample. It consists of the follwoing parts that are accompanied by a tutorial:

- [Package requirements](#package-requirements)
- [Parameter estimation](#parameter-estimation)
  This section explains how to do the parameter estimation from bulk WGS data. 
- [Simulated data](#simulated-data)
  This section explains how the simulated data were generated and how the model was learnt on them.
- [Published single cell WGS data](#published-single-cell-wgs-data)
  This section explains how the model was applied to published single-cell WGS data.
- [Real data](#real-data)
  This section explains how the model was applied to the actual data.

## Package requirements

Download and install R v4.2.0 and the following packages on your computer:
- CRAN
  - ape v5.6-2
  - phytools v1.2-0
  - phangorn v2.10.0
  - castor v1.7.5
  - TreeTools v1.8.0
  - deSolve v1.33
  - openxlsx v4.2.5
  - cdata v1.2.0

We have bundled R functions to model drift and selection in growing and homeostatic tissues in the R package FLORENCE. Please download and install from (https://github.com/VerenaK90/FLORENCE).

Download and install python3 and pyABC von your computer

## Simulated data

### Tree and VAF simulation

To test inference of drift and selection dynamics, we generated simulated division trees, generated simulated WGS data by binomial sampling of a simulated bulk and ran the model inference on these data. To reproduce the analysis, run the scripts [Neutral_tree.R](Simulated_data/Neutral_tree.R) and [Selected_tree.R](Simulated_data/Selected_tree.R) to simulate neutrally evolving trees and trees with a selected clone instigated at 20 years. The scripts source the scripts [Simulate_trees.R](Simulated_data/Simulate_trees.R) and [Tree_post_processing](Simulated_data/Tree_post_processing.R), which contain the modalities for tree simulation. Finally, to generate the VAF distributions which were used to assess the performance of FLORENCE, run [Study_design.R](Simulated_data/Study_design.R), which adds binomial noise to the simulated VAFs and generates the simulated data sets as outlined in (Metadata/Sample_information_simulated_data.xlsx).

### Parameter inference

To estimate the model parameters, we used FLORENCE in conjunction with pyABC. To re-run the analysis refer to the folder (Parameter_estimation) and modify Simulated_data/Run_model_sim*x.R according to the sample specification. Specifically, you should provide the following information:

- *patient.id*, the ID/name of the analyzed subject
- *age*, the age (in days)
- *snvs*, a named list containing data frames with VAF and Depth information for each individual 
- *depth*, the sequencing depth used to generate the data 
- *min.vaf*, the smallest VAF in the data that is to be compared to the model. Defaults to 0.05; we used 0.05 for all instances except for single cell WGS and simulated data with 270x, where we used 0.01.
- *min.clone.size*, the minimal clone size that can be detected by the model. Defaults to 0.05; we used 0.05 for all instances except for single cell WGS and simulated data with 270x, where we used 0.01.
- *min.prior.size*, the lower limit of clone sizes scanned by the parameter estimation. Parameter setzs associated with clones < min.clone.size will be evaluated with the neutral model. Defaults to 0.01; we used 0.01 for all instances except for single cell WGS and simulated data with 270x, where we used 0.001.
- *use.sensitivity* should sequencing sensitivity information be included in addition to binomial sampling? Defaults to F; if T, a matrix *false.negative.per.vaf* with columns corresponding to the measured VAFs and rows corresponding to individual measurements of the false negative rate at this VAF in addition to binomial noise must be provided. We used use.sensitivity=F.
  
### Analysis and plots

Upon parameter estimation, the posterior probabilities were analyzed using the script (Analyze_and_plot_fits_sim_data.R)[Analysis_and_figures/Analyze_and_plot_fits_sim_data.R]. This script 
- plots for each instance the model fit and the highest density estimates of each parameter
- computes statistics of the fits - %posterior probability supporting the selection model and the neutral model
- evaluates true and false positives for different clone sizes and sequencing depths, constructs the corresponding ROC curves and computes the AUC at the selected operating point
- this script generates the figure panels for **Fig. 1g-j** and **Fig. S1d-f**.


