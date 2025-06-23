# Central R script
# by Camille Sicangco
# Created 27 September 2024

# SETUP ########################################################################

# Load packages and functions
source("R/load_packages.R")
source("R/functions/data_processing_functions.R")
source("R/functions/analysis_functions.R")
source("R/functions/plotting_functions.R")

# Replace plantecophys::Photosyn with custom version to use the Heskel et al. 2017 R(T) equation
environment(Photosyn_custom) <- asNamespace("plantecophys")
assignInNamespace("Photosyn", Photosyn_custom, ns = "plantecophys")

# Replace calc_costgain to compute CGnet with netorig
environment(calc_costgain_netorig) <- asNamespace("gsthermal")
assignInNamespace("calc_costgain", calc_costgain_netorig, ns = "gsthermal")

# Replace calc_A with correction
environment(calc_A_corr) <- asNamespace("gsthermal")
assignInNamespace("calc_A", calc_A_corr, ns = "gsthermal")

# WTC DATA PROCESSING ##########################################################

# Process WTC4 input data
source("R/processing/WTC4_data_processing.R")

# Fit Tcrit, T50
source("R/processing/T50_fitting.R")

# Fit R-T response
source("R/processing/RT_fitting.R")

# Fit A-Ci T response
source("R/processing/ACi_T_fitting.R")

# ANALYSIS #####################################################################

## Theoretical simulations -----------------------------------------------------

# Force with increasing Tair and various Ps values
source("R/analysis/T_range_testing.R")

# Run instantaneous simulations to examine behaviour at specific Tair values
source("R/analysis/inst_sims.R")

## Sensitivity analysis --------------------------------------------------------

# Tcrit, T50
source("R/analysis/Tthreshold_sensitivity.R")

## WTC simulations -------------------------------------------------------------

# Fit gs models to WTC4 data
source("R/analysis/WTC_simulations.R")

## Supplementary figures -------------------------------------------------------
source("R/analysis/supporting_info.R")