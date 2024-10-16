This folder contains reproducible code for the simulations and real-data application in the paper:
"Sparse Causal Effect Estimation using Two-Sample Summary Statistics in the Presence of Unmeasured Confounding" by Huang et al. 2024.

# Simulations scripts:
- run_dgp1.R
- run_dgp2.R
- run_dgp3.R
These scripts save outputs in the `result` folder.

# Scripts to produce plots:
- summarize_dgp1.R
- summarize_dgp2.R
- summarize_dgp3.R
These scripts load the previous outputs in the `result` folder.

# Real-data file:
- data/expression.RData
See data/data_info.txt for data source information.

# Real-data application script:
- application_glp1r.R
