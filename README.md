# Clonal_hemopoiesis
This repository contains code to reproduce analysis in Körber et al., Quantifying drift and selection from a single bulk tissue sample. It consists of three parts that are acco:

- [Parameter estimation](#parameter-estimation)
  This section explains how to do the parameter estimation from bulk WGS data. 
- [Simulated data](#simulated-data)
  This section explains how the simulated data were generated and how the model was learnt on them.
- [Published single cell WGS data](#published-single-cell-wgs-data)
  This section explains how the model was applied to published single-cell WGS data.
- [Real data](#real-data)
  This section explains how the model was applied to the actual data.

## Simulated data

To test inference of drift and selection dynamics, we generated simulated division trees, generated simulated WGS data by binomial sampling of a simulated bulk and ran the model inference on these data. To reproduce the analysis, install the following packages:
- ape v5.6-2panied by tut
- phytools v1.2-0
- phangorn v2.10.0
- castor v1.7.5
- TreeTools v1.8.0
Run the scripts [Neutral_tree.R](Simulated_data/Neutral_tree.R) and [Selected_tree.R](Simulated_data/Selected_tree.R) to simulate neutrally evolving trees and trees with a selected clone instigated at 20 years. The scripts source the scripts [Simulate_trees.R](Simulated_data/Simulate_trees.R) and [Tree_post_processing](Simulated_data/Tree_post_processing.R), which contain the modalities for tree simulation. For a more detailed explanation on how to use these scripts, refer to Tutorial 1.

###


