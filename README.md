# FishCost


`FishCost` is a statistical framework that conducts a cost-effective analysis of commercial fisheries (fishery-dependent) and scientific survey (fishery-independent) data collection programs.


> Based on the paper [Getting more for less: benchmarking the cost-effectiveness of fishery-dependent and  -independent monitoring programs](https://www.google.com), currently under review in *Canadian Journal of Fisheries and Aquatic Sciences*.


Instructions
=============
Please make sure to install the following R packages in order to run `FishCost` 
* devtools
* TMB (Template Model Builder)
* gridConstruct

#### Installation 
    install.packages("devtools"); library("devtools")
    install.packages("TMB")
    devtools::install_github("kaskr/gridConstruct",subdir="gridConstruct")


Workflow
=========
The `FishCost` framework is built upon three consecutive R scripts:
* SimuAbu.R 
* AbuProcessing.R
* DEA.R

**Note: Scripts have to be run in the exact same order**

### >SimuAbu.R 
This script simulates fish abundance given a pre-defined sampling strategy. Abundance estimations are done through the LGNB spatio-temporal species distribution model developed by Rufener et al. (2021 - paper *in press*).

- Please refer to the [LGNB](https://github.com/mcruf/LGNB) GitHub repositroy for a full instruction on the LGNB model.


### >AbuProcessing.R
Once abundance were simulated, it is time to process the simulation results. For each sampling strategy and simulation therein, total abundances are calculated from the LGNB model output. The script is built in such way that it processes automatically the results according to the input data that is chosen by the user (survey, commercial, or integrated data).


### >DEA.R
This script conducts the Data Envelopment Analysis (DEA) based on data-specific sampling costs and the simulated abundances that were processed in the previous script.

