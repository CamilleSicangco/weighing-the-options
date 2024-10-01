# Central R script
# by Camille Sicangco
# Created 27 September 2024

# Load packages and functions
source("R/load_packages.R")

# Force models with increasing air temperature
source("R/analysis/T_range_testing.R")

# Process WTC4 data
source("R/processing/WTC4_data_processing.R")

# Fit gs models to WTC4 data
source("R/analysis/run_model_WTC4.R")