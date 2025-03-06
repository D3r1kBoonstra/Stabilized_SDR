# Stabilized SDR (SSDR)

## Overview

This repository accompanies the research paper:

**"Precision Matrix Regularization in Sufficient Dimension Reduction for Improved Quadratic Discriminant Classification"** 

*Authors:* Derik T. Boonstra, Rakheon Kim, and Dean M. Young

Stabilized SDR (SSDR), is a multiclass and distribution-free Sufficient Dimension Reduction (SDR) method that employs user-specified precision-matrix shrinkage estimators to stabilize the projection-matrix and supervised classifier. 

## Repository Structure

- `datasets/`: Contains datasets used for real data applications.
- `mc_sims/`: R Scripts for Monte Carlo simulations.
- `real_data_applications/`: Scripts for real data applications. 
- `saved_sims/`: Saved results from simulations and applications.
- `Requirements.R`: R script listing required packages. Run before any computations. 
- `Stabilized_SDR.Rproj`: RStudio project file. To efficiently run code, open project first. 
- `figures.R`: Script for generating figures from the paper.
- `wrapper_fns.R`: Wrapper functions used in the analysis.
  
Note: Analysis is primarily in R. Matlab is used for the ENDS methods only.

If you find this work useful, please cite our paper:

@article{boonstra2025,  
  author = {Boonstra, Derik T. and Kim, Rakheon and Young, Dean M.},  
  title = {Precision Matrix Regularization in Sufficient Dimension Reduction for Improved Quadratic Discriminant Classification},  
  journal = {tbd},  
  year = {tbd},    
}
