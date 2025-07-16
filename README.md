# Weighing the options
Camille K. Sicangco, Manon E. B. Sabot, John E. Drake, Mark G. Tjoelker, & Belinda E. Medlyn

## Overview
Repository containing the code to reproduce the results in Sicangco et al., "Weighing the options: a test of alternative stomatal optimisation models at high temperatures". Submitted to New Phytologist.

## Instructions
To use this code, you will need to install `gsthermal`, a custom R package built for this analysis. You can find the package here: https://github.com/CamilleSicangco/gsthermal.

To regenerate our results, run `R/central_script.R`. This script will guide you through loading the required packages and functions, processing the data, and conducting the analysis. 

## WTC dataset
The models are evaluated against the dataset from Drake et al. (2018), "Trees tolerate an extreme heatwave via sustained transpirational cooling and increased leaf thermal tolerance". The data used for our analysis can be accessed freely at:

* Fluxes, leaf temperatures, water potentials: [Drake et al. (2018)](http://doi.org/10.4225/35/5a36f61f150f3)
* A-Ci T-response parameters: [PPC-TGlob_V1.0, Kumarathunge et al. (2019)](https://bitbucket.org/Kumarathunge/aci-tglob/src/master/)
* Soil water content, raw Fv/Fm, R-T: [data/in/raw](https://github.com/CamilleSicangco/ThermalCost/tree/main/data/in)

## Contact
Camille Sicangco: C.Sicangco@westernsydney.edu.au