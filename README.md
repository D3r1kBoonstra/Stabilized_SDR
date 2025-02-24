# SDR Shrinkage (SDRS)

## Overview

This repository accompanies the research paper:

**"Sufficient Dimension Reduction Shrinkage for Improved Quadratic
Discriminant Classification"** 

*Authors:* Derik T. Boonstra, Rakheon Kim, Gabriel J. Odom, and Dean M. Young

SDR Shrinkage (SDRS), is a multiclass distribution-free Sufficient Dimension Reduction (SDR) method that employs user-specified precision-matrix shrinkage estimators to stabilize both the projection-matrix and supervised classifier. 

## Repository Structure

- `datasets/`: Contains datasets used for real data applications.
- `mc_sims/`: R Scripts for Monte Carlo simulations.
- `real_data_applications/`: Scripts for real data applications. 
- `saved_sims/`: Saved results from simulations and applications.
- `Requirements.R`: R script listing required packages. Run before any computations. 
- `SDRshrinkage.Rproj`: RStudio project file. To efficiently run code, open project first. 
- `figures.R`: Script for generating figures from the paper.
- `wrapper_fns.R`: Wrapper functions used in the analysis.
  
Note: Analysis is primarily in R. Matlab is used for the ENDS methods only.

If you find this work useful, please cite our paper:

@article{boonstra2025,  
  author = {Boonstra, Derik T. and Kim, Rakheon and Odom, Gabriel J. and Young, Dean M.},  
  title = {Sufficient Dimension Reduction Shrinkage for Improved Quadratic Discriminant Classification},  
  journal = {tbd},  
  year = {tbd},    
}
