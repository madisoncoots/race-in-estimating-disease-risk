# Author: Madison Coots
# Date: November 5, 2023
# ======================
# Code for all disease risk predictions. 
# Coots et al. (2023)

library(tidyverse)
library(rstudioapi)
library(here)

directory_path <- dirname(getActiveDocumentContext()$path)
directory_path %>% str_c("/diabetes/diabetes_risk_models.R")

diabetes_script_path <- directory_path %>% str_c("/diabetes/diabetes_risk_models.R")
cvd_script_path <- here(directory_path, "cvd/2013_pce.R")
bc_script_path <- here(directory_path, "bc/gail_equation.R")

setwd(directory_path)
# =============================================================================
# Run risk prediction script for each disease to generate and save data objects
# with risk predictions to the data/processed/ directory of this repository.
# =============================================================================
source(diabetes_script_path)

source(cvd_script_path)

source(bc_script_path)
