# Central R script
# by Camille Sicangco
# Created 27 September 2024

# Load packages and functions
source("R/load_packages.R")
source("R/functions/data_processing_functions.R")
source("R/functions/analysis_functions.R")

# Replace plantecophys::Photosyn with custom version 
# to use the Heskel et al. 2017 R(T) equation
environment(Photosyn_custom) <- asNamespace("plantecophys")
environment(Photosyn_custom2) <- asNamespace("plantecophys")
assignInNamespace("Photosyn", Photosyn_custom, ns = "plantecophys")

# Replace calc_costgain to compute CGnet with netorig
environment(calc_costgain_netorig) <- asNamespace("gsthermal")
assignInNamespace("calc_costgain", calc_costgain_netorig(), ns = "gsthermal")


# Force models with increasing air temperature
source("R/analysis/T_range_testing.R")

# Process WTC4 data
source("R/processing/WTC4_data_processing.R")

# Fit gs models to WTC4 data
source("R/analysis/run_model_WTC4.R")