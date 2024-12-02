# Author: Madison Coots
# Date: October 31, 2023
# ----------------------
# This file makes breast cancer risk predictions.
# Coots et al. (2023)

library(readr)
library(BCRA) # Library for breast cancer risk
library(tidyverse)
library(nnet) # for multinomial regression
library(here)
library(rstudioapi)

working_directory <- getwd()
if (basename(working_directory) == "code") {
  # being called from master script
  directory_path <- here(working_directory, 'bc/')
} else {
  # being run on its own
  directory_path <- dirname(getActiveDocumentContext()$path)
}

data_object_read_path <- directory_path %>% str_remove("code/bc") %>% str_c("data/bc/processed/bc_data.rds")
data_object_write_path <- directory_path %>% str_remove("code/bc") %>% str_c("data/bc/processed/")

data <- readRDS(data_object_read_path)

# ======================================================================
# Predicting race-aware 5-year breast cancer risk (standard)
# ======================================================================

all_data_with_race_aware_risk <- data %>%
  mutate(race_aware_risk = absolute.risk(data))

# ======================================================================
# Predicting race-blind 5-year breast cancer risk (non-standard)
# ======================================================================

# ====== Need to train multinomial model to compute P(R = r | X) =======
race_weight_model_formula <- Race ~ T1 + T2 + US_born + AgeMen +
  Age1st + N_Biop + HypPlas + N_Rels

race_weight_model <- multinom(race_weight_model_formula, 
                              data = data, 
                              weights = wtmec8yr)

race_weights <- predict(race_weight_model, data, type = "probs") %>%
  data.frame() %>% 
  rename(white_weight = X1,
         black_weight = X2,
         us_hispanic_weight = X3,
         non_us_hispanic_weight = X5,
         asian_weight = X11) %>%
  mutate(ID = data$ID)

# ======================================================================

all_data_with_race_blind_risk <- data %>%
  mutate(white_risk = absolute.risk(data %>% mutate(Race = 1)),
         black_risk = absolute.risk(data %>% mutate(Race = 2)),
         us_hispanic_risk = absolute.risk(data %>% mutate(Race = 3)),
         non_us_hispanic_risk = absolute.risk(data %>% mutate(Race = 5)),
         asian_risk = absolute.risk(data %>% mutate(Race = 11))) %>%
  left_join(race_weights, by = c("ID")) %>%
  mutate(race_blind_risk =
           white_weight * white_risk +
           black_weight * black_risk +
           us_hispanic_weight * us_hispanic_risk +
           non_us_hispanic_weight * non_us_hispanic_risk +
           asian_weight * asian_risk) %>%
  select(-white_risk, -black_risk, -us_hispanic_risk, -non_us_hispanic_risk, -asian_risk,
         -white_weight, -black_weight, -us_hispanic_weight, -non_us_hispanic_weight, -asian_weight)

# ======================================================================
# Saving data and both predictions
# ======================================================================

all_data_with_both_risk <-
  all_data_with_race_aware_risk %>%
  left_join(all_data_with_race_blind_risk %>%
              select(ID, race_blind_risk), by = c("ID")) %>%
  mutate(race_aware_risk = race_aware_risk / 100,
         race_blind_risk = race_blind_risk / 100)

saveRDS(all_data_with_both_risk, file = paste(data_object_write_path, "bc_data_with_risk.rds", sep = ""))
