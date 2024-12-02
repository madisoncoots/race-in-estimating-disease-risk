# Author: Madison Coots
# Date: November 5, 2023
# ======================
# Code for all data processing. 
# Coots et al. (2023)

library(tidyverse)
library(rstudioapi)
library(here)

directory_path <- dirname(getActiveDocumentContext()$path)

diabetes_script_path <- here(directory_path, "diabetes/make_data_diabetes.R")
cvd_script_path <- here(directory_path, "cvd/make_data_cvd.R")
bc_script_path <- here(directory_path, "bc/make_data_bc.R")

setwd(directory_path)
# =============================================================================
# Run make_data_*.R script for each disease to generate and save data objects
# to the data/ directory of this repository
# =============================================================================

source(diabetes_script_path)

source(cvd_script_path)

source(bc_script_path)
