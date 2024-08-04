# Clonal_hematopoiesis
This repository contains code to reproduce analysis in Körber et al., Detecting and quantifying clonal selection in somatic stem cells. It consists of the following parts:

- [Before you start](#before-you-start)
  This section gives an overview of the software and data used for the analysis.
- [Simulated data](#simulated-data)
  This section explains how the simulated data were generated and how the model was applied to these data.
- [Published single cell WGS data](#published-single-cell-wgs-data)
  This section explains how the model was applied to published single-cell WGS data.
- [Bulk WGS data from human heme](#bulk-WGS-data-from-human-heme)
  This section explains how the model was applied to the newly generated bulk WGS data from human bone marrow.
- [Bulk WGS data from human brain](#bulk-WGS-data-from-human-brain)
  This section explains how the model was applied to published bulk WGS data from human brain.
  
To reproduce (parts of) the analysis, please download the supplementary tables, the associated data from Mendeley (doi: 10.17632/yvxdb7t3yk.1) and adjust the directories in [Settings.R](Settings.R). Refer to [Before you start](#before-you-start) for a list of R libraries that should be installed.

## Before you start

Download and install R v4.2.0 and the following packages on your computer:
- CRAN
  - bedr v1.0.7
  - ape v5.6-2
  - phytools v1.2-0
  - phangorn v2.10.0
  - castor v1.7.5
  - TreeTools v1.8.0
  - deSolve v1.33
  - openxlsx v4.2.5
  - cdata v1.2.0
  - ggpubr v0.4.0
  - phangorn v2.10.0
  - RRphylo v2.7.0
  - ggplot2 v3.4.2
  - cgwtools v3.3
  - ggVennDiagram v1.2.2
  - ggbeeswarm v0.6.0
  - ggformula v0.10.2
  - ggsci v2.9
  - Hmisc v4.7.1
  - lemon v0.4.5
  - data.table v1.14.2
  - RColorBrewer v1.1.3
  - ggridges v0.5.4
  - wesanderson v0.3.6
  - HDInterval v0.2.2
  - reshape2 v1.4.4
  - dplyr v1.0.9
  - scales v1.2.1

We have bundled R functions to model drift and selection in growing and homeostatic tissues in the R package SCIFER. Please download and install from (https://github.com/VerenaK90/SCIFER/tree/paper) and refer to the accompanying vignette for further information on SCIFER.

Download and install python3 and pyABC von your computer.

Associated data can be found on Mendeley (doi: 10.17632/yvxdb7t3yk.1).

To reproduce the analysis, please check the file structure and make sure the required scripts and input files are in the correct location.

## Variant allele frequencies shaped by drift and selection in homeostatic tissues

SCIFER predicts that drift and selection shape the variant allele frequency distribution distinctly. We predict VAF distributions for different parameter sets and time points in the script [Theoretical_model_performance.R](Analysis_and_figures/Theoretical_model_performance.R). The script also produces the figure panels in **Figure 1** and **Supplementary Fig. 1**.

## Simulated data

### Tree and VAF simulation

To test the ability of SCIFER to distinguish selection and drift dynamics, we simulated division trees. Of these, we simulated WGS data by binomial sampling of a simulated bulk sample and applied SCIFER to these data. To reproduce the analysis, run the scripts [Neutral_tree.R](Simulated_data/Neutral_tree.R) and [Selected_tree.R](Simulated_data/Selected_tree.R) to simulate neutrally evolving trees and trees with a selected clone instigated at 20 years. To generate the VAF distributions which were used to assess the performance of SCIFER, run [Study_design.R](Simulated_data/Study_design.R), which adds binomial noise to the simulated VAFs and generates the simulated data sets as outlined in (MetaData/Sample_information_simulated_data.xlsx). To simulate a phylogenetic tree with both stem and progenitor cells as described in Supplementary Note 1, run the script [Neutral_tree_S_P.R](Simulated_data/Neutral_tree_S_P.R). This script also produces the plots shown in Supplementary Figure 1.

### Parameter estimation

To estimate the model parameters, we used SCIFER in conjunction with pyABC. To re-run the analysis refer to the folder (Parameter_estimation) and modify Parameter_estimation/Run_model_sim*x.R according to the sample specification. Specifically, you should provide the following information:

- `patient.id`, the ID/name of the analyzed subject
- `age`, the age (in days, can be retrieved from Sample_information_simulated_data.xlsx)
- `snvs`, a named list containing data frames with VAF and Depth information for each individual (e.g., as provided in RData/Simulated_data/SNVs_90x.RData)
- `depth`, the sequencing depth used to generate the data 
- `min.vaf`, the smallest VAF in the data that is to be compared to the model. Defaults to 0.05; we used `min.vaf=0.05` for all instances except for single cell WGS and simulated data with 270x, where we used `min.vaf=0.01`.
- `min.clone.size`, the minimal clone size that can be detected by the model. Defaults to 0.05; we used `min.clone.size=0.05` for all instances except for single cell WGS and simulated data with 270x, where we used `min.clone.size=0.01`.
- `min.prior.size`, the lower limit of clone sizes scanned by the parameter estimation. Parameter sets associated with clones < min.clone.size will be evaluated with the neutral model. Defaults to 0.01; we used `min.prior.size=0.01` for all instances except for single cell WGS and simulated data with 270x, where we used `min.prior.size=0.001`.
- `use.sensitivity` should sequencing sensitivity information be included in addition to binomial sampling? Defaults to `F`; if `T`, a matrix `false.negative.per.vaf` with columns corresponding to the measured VAFs and rows corresponding to individual measurements of the false negative rate at this VAF in addition to binomial noise must be provided. We used `use.sensitivity=F`.

The script Run_model_sim*x.R is to be sourced by [ABC_fit.py](Parameter_estimation/ABC_fit.py). Hence, please also adjust the input and output paths in this file. The python script also contains the definition of the prior distributions. Note that we chose priors running between 0 and 0.99 for `s` and between 0 and 1 for `t_s`, but these values are relative values only that will be converted into absolute values in the model script [Bayesian_fit.R](Parameter_estimation/Bayesian_fit.R). Specifically, the minimal and maximal values of `s` and `t_s` are chosen such that the clone has a minimal size of `min.prior.size` and a maximal size of 1. Analogous conversion of these two parameters into absolute values is also necessary when inspecting the parameter estimates later on (refer to [Analyze_and_plot_fits_sim_data.R](Analysis_and_figures/Analyze_and_plot_fits_sim_data.R)).
  
### Analysis and plots

Upon parameter estimation, the posterior probabilities (which can be downloaded from Mendeley) were analyzed using the script [Analyze_and_plot_fits_sim_data.R](Analysis_and_figures/Analyze_and_plot_fits_sim_data.R). This script 
- plots for each instance the model fit and the highest density estimates of each parameter
- computes statistics of the fits - %posterior probability supporting the selection model and the neutral model
- evaluates true and false positives for different clone sizes and sequencing depths, constructs the corresponding ROC curves and computes the AUC at the selected operating point
- this script generates the figure panels for **Fig. 2**.

## Published single cell WGS data

We tested our model on pseudo-bulk WGS data from published single-cell WGS data of three studies, Lee-Six et al., Nature, 2018, Mitchell et al., Nature, 2022 and Fabre et al., 2022. In a first step, we generated the pseudo-bulks for each study as outlined in the following.

### Generation of pseudo-bulk data 

#### Lee-Six et al.
Download the data from https://doi.org/10.17632/yzjw2stk7f.1 and store them in ./Published_data/Lee-Six_et_al/Caveman/. In addition, download the re-called data from our repository (https://doi.org/10.17632/yvxdb7t3yk.1) and store them in ./Lee-Six_et_al/Mutect_Strelka.

The script [Pseudo_VAFs_LeeSix_et_al.R](Data_preprocessing/Pseudo_VAFs_LeeSix_et_al.R) first compares the results of Caveman and Mutect/Strelka (**Extended Data Fig. 3b,c**). Thereafter, it generates the pseudo-bulk data and plots the VAF distributions as shown in **Fig. 3b, Extended Data Fig. 3d**. Finally, it stores two list objects containing the VAFs in [Caveman/SNVs.RData](RData/Lee-Six_et_al/Caveman/SNVs.RData) and [Mutect_Strelka/SNVs.RData](RData/Lee-Six_et_al/Mutect_Strelka/SNVs.RData). Alternatively, these objects can be directly downloaded from Mendeley.

#### Mitchell et al.

Download the data from https://data.mendeley.com/datasets/np54zjkvxr/1 and structure them like: ./Published_data/Mitchell_et_al/*/.

The script [Pseudo_VAFs_Mitchell_et_al.R](Data_preprocessing/Pseudo_VAFs_Mitchell_et_al.R) generates the pseudo-bulk data and plots the VAF distributions as shown in **Fig. 3h, n, Extended Data Fig. 3e, f**. It also plots the trees as shown in **Fig. 3g, m**. Finally, it stores a list object containing the VAFs for each sample in (SNVs.RData)[RData/Mitchell_et_al/SNVs.RData]. Alternatively, this object can be directly downloaded from Mendeley.

#### Fabre et al.
Download the data of id2259 from doi.org/10.6084/m9.figshare.15029118 and structure them like: ./Published_data/Fabre_et_al/*/.

The script [Pseudo_VAFs_Fabre_et_al.R](Data_preprocessing/Pseudo_VAFs_Fabre_et_al.R) generates the pseudo-bulk data and plots the VAF distribution as shown in **Extended Data Fig. 3g**. It also plots the trees as shown in **Extended Data Fig. 3g**. Finally, it stores a list object containing the VAFs in [SNVs.RData](RData/Fabre_et_al/SNVs.RData). Alternatively, this object can be directly downloaded from Mendeley.

### Parameter estimation

As with the simulated data, we used SCIFER in conjunction with pyABC to estimate parameters. To re-run the analysis refer to the folder [Parameter_estimation](Parameter_estimation) and modify [Run_model_scWGS.R](Run_model_scWGS.R) according to the sample specification, with emphasis on the following information

- `patient.id`, the ID/name of the analyzed subject
- `age`, the age (in days; can be retrieved from MetaData/Sample_info_published_data.xlsx)
- `snvs`, a named list containing VAF information for each individual (e.g., provided in RData/Mitchell_et_al/SNVs.RData)
- `depth`, the sequencing depth used to generate the data; can be arbitrarily set here, as we will use the "sc" mode (see below, `seq.type`)
- `min.vaf`, the smallest VAF in the data that is to be compared to the model. Defaults to 0.05; we used `min.vaf=0.01` due to the high pseudo-bulk coverage.
- `min.clone.size`, the minimal clone size that can be detected by the model. Defaults to 0.05; we used `min.clone.size=0.01` due to the high pseudo-bulk coverage.
- `min.prior.size`, the lower limit of clone sizes scanned by the parameter estimation. Parameter setzs associated with clones < min.clone.size will be evaluated with the neutral model. We used `min.prior.size=0.001` due to the high pseudo-bulk coverage.
- `ncells`, the number of sequenced cells
- `seq.type`, has to be set to `seq.type="sc"`, as we analyze pseudo-bulks from single-cells
- `use.sensitivity` should sequencing sensitivity information be included in addition to binomial sampling? Defaults to F; if T, a matrix `false.negative.per.vaf` with columns corresponding to the measured VAFs and rows corresponding to individual measurements of the false negative rate at this VAF in addition to binomial noise must be provided. We set `use.sensitivity=F`.

In comparison to bulk WGS data sequenced at 90x, we here lowered the lower VAF limits, as the single-cell sequencing data provide higher resolution. Moreover, we specified the number of sequenced cells and run the parameter estimation in single-cell mode, allowing for the simulation of single-cell sequencing and generation of pseudo-bulk data thereof. 

The script [Run_model_scWGS.R](Parameter_estimation/Run_model_scWGS.R) is to be sourced by [ABC_fit.py](Parameter_estimation/ABC_fit.py). Hence, please also adjust the input and output paths in this file. The python script also contains the definition of the prior distributions. Note that we chose priors running between 0 and 0.99 for s and between 0 and 1 for t_s, but these values are relative values only that will be converted into absolute values by the model script [Bayesian_fit.R](Parameter_estimation/Bayesian_fit.R). Specifically, the minimal and maximal values of s and t_s are chosen such that the clone has a minimal size of `min.prior.size and a maximal size of 1. Analogous conversion of these two parameters into absolute values is also necessary when inspecting the parameter estimates later on (refer to [Plot_fits_published_data.R](Analysis_and_figures/Plot_fits_published_data.R)).

#### Parameter estimation for a 2-clone model for sample KX004

In a second step, we used SCIFER in conjunction with pyABC to estimate parameters for a two-clone model for sample KX004, where multiple expanded clades were detected in the phylogenetic tree. Based on the tree phylogeny, the model fit was performed assuming that the two clones originated via `branched` evolution. To re-run the analysis refer to the folder [Parameter_estimation](Parameter_estimation) and modify [Run_model_scWGS_2_branched_clones.R](Simulated_data/Run_model_scWGS_data_2_branched_clones.R).

The script Run_model_scWGS_2_branched_clones.R is to be sourced by [ABC_fit_branched.py](Parameter_estimation/ABC_fit_branched.py). Hence, please adjust the input and output paths in this file. The python script also contains the definition of the prior distributions. Note that for the 2-clone model, we parametrized the model with prior distributions on the 2 relative size variables `size1` and `size2`, informed by the results obtained with the one-clone model (see Supplementary Table 10 for details). The relative values will be converted into absolute values for `s1` and `s2` in the model script [Bayesian_fit_multiclone.R](Parameter_estimation/Bayesian_fit_multiclone.R) (see Methods for details). Analogous conversion of these two parameters into absolute values is also necessary when inspecting the parameter estimates later on (refer to [Plot_fits_KX004_SCIFER2.R](Analysis_and_figures/Plot_fits_KX004_SCIFER2.R)). 

### Analysis and plots

Once the parameter estimation has been finished, extract .csv-files from the .db files using the function abc-export (or directly download them from Mendeley). The fits can then be inspected using the script [Plot_fits_published_data.R](Analysis_and_plots/Plot_fits_published_data.R).

This script
- plots model vs data for each sample (as shown in **Fig. 3b,h,n, Extended Data Figs. 3d,f,g**)
- classifies individual cases as neutrally evolving or selected (as shown in **Fig. 3c,i,o**)
- computes highest density estimates for the parameters (as shown in **Figs. 3d-f,k,l,p**)
- compares the estimated age of the selected clones to the original studies (as shown in **Fig. 3j**).

#### Analyzing the results of a two-clone model for sample KX004

Upon parameter estimation, estract a .csv-file from the .db-file using the funciton abc-export (or directly download the file from Mendeley). The fit can then be inspected using the script [Plot_fits_KX004_SCIFER2.R](Analysis_and_figures/Plot_fits_KX004_SCIFER2.R), which 
- plots model vs data for sample KX004 (as shown in **Extended Data Fig. 3j**)
- plots the posterior probability for the first and second selected clone obtained with the 2-clone model (as shown in **Extended Data Fig. 3j**)
- computes highest density estimates for the parameters (as shown in **Extended Data Fig. 3k**)
- compares the estimated parameters with the results from the one-clone model (**Extended Data Fig. 3k**)

## Bulk WGS data from human heme

### Data pre-processing

We ran the model on the filtered SNVs (called using Strelka and Mutect2, see manuscript for details), which we stored, for convenience, in the list object (./RData/WGS_heme/SNVs.RData). 

### Parameter estimation 
#### One-clone model

As with the other data types, we used SCIFER in conjunction with pyABC to estimate parameters. To re-run the analysis refer to the folder [Parameter_estimation](Parameter_estimation) and modify [Run_model_heme_WGS_data.R](Simulated_data/Run_model_heme_WGS_data.R) according to the sample specification, with emphasis on the following information

- `patient.id`, the ID/name of the analyzed subject
- `sort`, the cell sort to be analyzed ("CD34", "CD34_deep", "MNC", "MNC_minus_T" or "PB_gran")
- `age`, the age (in days)
- `snvs`, a named list containing a data frame with VAFs and depths for each individual, as provided in RData/WGS_heme/SNVs.RData 
- `depth`, the sequencing depth used to generate the data 
- `min.vaf`, the smallest VAF in the data that is to be compared to the model. We used `min.vaf=0.05` for 90x/120x WGS and `min.vaf=0.02`for 270x WGS, according to the respective detection limits.
- `min.clone.size`, the minimal clone size that can be detected by the model. We used `min.clone.size=0.05` for 90x/120x WGS and `min.vaf=0.02` for 270x WGS, according to the respective detection limits.
- `min.prior.size`, the lower limit of clone sizes scanned by the parameter estimation. Parameter setzs associated with clones < min.clone.size will be evaluated with the neutral model. We used `min.prior.size=0.01` for 90x/120x WGS and `min.prior.size=0.001` for 270x WGS, according to the respective detection limits.
- `use.sensitivity` should sequencing sensitivity information be included in addition to binomial sampling? Defaults to F; if T, a matrix `false.negative.per.vaf` with columns corresponding to the measured VAFs and rows corresponding to individual measurements of the false negative rate at this VAF in addition to binomial noise must be provided. We used `use.sensitivity=F`.

The script Run_model_heme_WGS.R is to be sourced by [ABC_fit.py](Parameter_estimation/ABC_fit.py). Hence, please also adjust the input and output paths in this file. The python script also contains the definition of the prior distributions. Note that we chose priors running between 0 and 0.99 for s and between 0 and 1 for t_s, but these values are relative values only that will be converted into absolute values by the model script [Bayesian_fit.R](Parameter_estimation/Bayesian_fit.R). Specifically, the minimal and maximal values of s and t_s are chosen such that the clone has a minimal size of `min.prior.size` and a maximal size of 1. Analogous conversion of these two parameters into absolute values is also necessary when inspecting the parameter estimates later on (refer to [Assess_fits_heme_WGS_data.R](Analysis_and_figures/Assess_fits_heme_WGS_data.R)).

To re-run the analysis for a model modification without size compensation for the selected clone (i.e., the overall stem cell pool expands according to the expansion of the selected clone) refer to the folder [Parameter_estimation](Parameter_estimation) and modify [Run_model_heme_WGS_data_nsc.R](Simulated_data/Run_model_heme_WGS_data_nsc.R) instead of [Run_model_heme_WGS_data.R](Simulated_data/Run_model_heme_WGS_data.R) and perform the model fits in analogy.

#### Two-clone model

In a second step, we used SCIFER in conjunction with pyABC to estimate parameters for a two-clone model. Here, we only used 270x WGS data from sorted CD34+ cells for parameter estimation. For each sample, two model fits were performed, assuming that the two clones originated via `branched` evolution vs `linear`evolution. To re-run the analysis refer to the folder [Parameter_estimation](Parameter_estimation) and modify [Run_model_heme_WGS_data_2_branched_clones.R](Simulated_data/Run_model_heme_WGS_data_2_branched_clones.R) and [Run_model_heme_WGS_data_2_linear_clones.R](Simulated_data/Run_model_heme_WGS_data_2_linear_clones.R) according to the sample specification, with emphasis on the following information

- `patient.id`, the ID/name of the analyzed subject
- `age`, the age (in days)
- `snvs`, a named list containing a data frame with VAFs and depths for each individual, as provided in RData/WGS_heme/SNVs.RData 
- `depth`, the sequencing depth used to generate the data 
- `min.vaf`, the smallest VAF in the data that is to be compared to the model. We used `min.vaf=0.02`, according to the detection limit of 270x WGS.
- `min.clone.size`, the minimal clone size that can be detected by the model. We used `min.clone.size=0.01`, according to the detection limit of 270x WGS.
- `min.prior.size`, the lower limit of clone sizes scanned by the parameter estimation. Parameter setzs associated with clones < min.clone.size will be evaluated with the neutral model. We used `min.prior.size=0.001`, according to the detection limit of 270x WGS.
- `use.sensitivity` should sequencing sensitivity information be included in addition to binomial sampling? Defaults to F; if T, a matrix `false.negative.per.vaf` with columns corresponding to the measured VAFs and rows corresponding to individual measurements of the false negative rate at this VAF in addition to binomial noise must be provided. We used `use.sensitivity=F`.

The scripts Run_model_heme_WGS_data_2_branched_clones.R and Run_model_heme_WGS_data_2_linear_clones.R are to be sourced by [ABC_fit_branched.py](Parameter_estimation/ABC_fit_branched.py) and [ABC_fit_linear.py](Parameter_estimation/ABC_fit_linear.py), respectively. Hence, please also adjust the input and output paths in these files. The python scripts also contains the definition of the prior distributions. Note that for the 2-clone model, we parametrized the model with prior distributions on the 2 relative size variables `size1` and `size2`, informed by the results obtained with the one-clone model (see Supplementary Table 10 for details). The relative values will be converted into absolute values for `s1` and `s2` in the model script [Bayesian_fit_multiclone.R](Parameter_estimation/Bayesian_fit_multiclone.R) (see Methods for details). Analogous conversion of these two parameters into absolute values is also necessary when inspecting the parameter estimates later on (refer to [Assess_fits_heme_WGS_data_2_clone_model.R](Analysis_and_figures/Assess_fits_heme_WGS_2_clone_model.R)).

### Analysis and plots

As with the pseudo-bulk data, extract .csv-files from the .db files using the function abc-export (or directly download them from Mendeley). The fits can then be inspected using the script [Plot_parameters_heme_WGS.R](Analysis_and_figures/Plot_parameters_heme_WGS.R).

This script sources the files [Assess_fits_heme_WGS_data.R](Analysis_and_figures/Assess_fits_heme_WGS_data.R), [Assess_fits_heme_WGS_data_2_clone_model.R](Analysis_and_figures/Assess_fits_heme_WGS_2_clone_model.R) and [Assess_fits_heme_WGS_data_no_size_compensation.R](Analysis_and_figures/Assess_fits_heme_WGS_data_no_size_compensation.R) to integrate/compare model fits between different models. Moreover, it
- plots data only for the neutral cases (**Fig. 5a**)
- plots model vs data for each sample (as shown in **Figs. 4a/b, 5d, Extended Data Fig. 6, Extended Data Fig. 7a-c,e, Extended Data Fig. 8a,b**)
- classifies individual cases as neutrally evolving or selected (as shown in **Figs. 4c-e, 5b,c, Supplementary Fig. 3**)
- compares estimated clone sizes to VAFs of known drivers (**Fig. 4f**)
- computes highest density estimates for the parameters (as shown in **Fig. 4g/h, 5e-h, 7a-c, Extended Data Fig. 7e, Fig.6a-c, Extended Data Fig. 9**)
- compares the estimated parameters between samples with and without selection (as shown in **Fig. 5i-k**)
- computes the cumulative age-incidence of CH driver acquisition (**Fig. 6d**)

## Bulk WGS data from human brain

This section describes the analysis of published bulk WGS data from Bae et al., Science, 2022 with SCIFER. 

### Data pre-processing

Download the supplementary tables (S1, S2, S6 and S7) from https://doi.org/10.1126/science.abm6222 and store them into the folder "Published_data/Bae_et_al". For our analysis, we used cases with a coverage >100x and, in addition, all samples with coverage >30x from case N7 (with multiple known driver mutations) as a positive control. Details on the samples included in our analysis can be found in Supplementary Table 8.

The script [Preprocess_bae_et_al.R](Data_preprocessing/Preprocess_bae_et_al.R) reads in the published data from Bae et al., and stores a list object, containing the SNVs, including VAF and depth, identified in each sample in the original publication in (SNVs_brain.RData)[RData/Bae_et_al/SNVs_brain.RData]. Alternatively, this object can be directly downloaded from Mendeley.

### Parameter estimation 

As with heme WGS data, we used SCIFER in conjunction with pyABC to estimate parameters. To re-run the analysis refer to the folder [Parameter_estimation](Parameter_estimation) and modify [Run_model_brain_WGS_data.R](Simulated_data/Run_model_brain_WGS_data.R) according to the sample specification, with emphasis on the following information

- `patient.id`, the ID/name of the analyzed subject
- `age`, the age (in days)
- `snvs`, a named list containing a data frame with VAFs and depths for each individual, as provided in RData/WGS_heme/SNVs.RData 
- `depth`, the sequencing depth used to generate the data 
- `min.vaf`, the smallest VAF in the data that is to be compared to the model. We used `min.vaf=0.05` for <150x and `min.vaf=0.02` for ≥150x WGS data, according to the respective detection limits.
- `min.clone.size`, the minimal clone size that can be detected by the model. We used `min.clone.size=0.05` for <150x and `min.vaf=0.02` for ≥150x WGS data, according to the respective detection limits.
- `min.prior.size`, the lower limit of clone sizes scanned by the parameter estimation. Parameter setzs associated with clones < min.clone.size will be evaluated with the neutral model. We used `min.prior.size=0.01` for <150x and `min.prior.size=0.001` for ≥150x WGS data, according to the respective detection limits.
- `use.sensitivity` should sequencing sensitivity information be included in addition to binomial sampling? Defaults to F; if T, a matrix `false.negative.per.vaf` with columns corresponding to the measured VAFs and rows corresponding to individual measurements of the false negative rate at this VAF in addition to binomial noise must be provided. We used `use.sensitivity=F`.

The script Run_model_heme_WGS.R is to be sourced by [ABC_fit.py](Parameter_estimation/ABC_fit.py). Hence, please also adjust the input and output paths in this file. The python script also contains the definition of the prior distributions. Note that we chose priors running between 0 and 0.99 for s and between 0 and 1 for t_s, but these values are relative values only that will be converted into absolute values by the model script [Bayesian_fit.R](Parameter_estimation/Bayesian_fit.R). Specifically, the minimal and maximal values of s and t_s are chosen such that the clone has a minimal size of `min.prior.size` and a maximal size of 1. Analogous conversion of these two parameters into absolute values is also necessary when inspecting the parameter estimates later on (refer to [Plot_fits_brain_WGS_data.R](Analysis_and_figures/Plot_fits_brain_WGS.R)).

### Analysis and plots

As with the heme data, extract .csv-files from the .db files using the function abc-export (or directly download them from Mendeley). The fits can then be inspected using the script [Plot_fits_brain_WGS_data.R](Analysis_and_figures/Plot_fits_brain_WGS_data.R).

This script 
- plots model vs data for each sample (as shown in **Figs. 8a-c**)
- classifies individual cases as neutrally evolving or selected (as shown in **Fig. 8d**)
- compares the number of cases with evidence for clonal selection across phenotype, sex and brain region (**Extended Data Fig. 10b**)
- computes the incidence of selection across age groups (**Fig. 8e**)
- computes highest density estimates for the parameters (as shown in **Fig. 8f,g, Extended Data Fig. 10c,d**)
- computes the cumulative age-incidence of the onset of clonal selection (**Fig. 8h**)

