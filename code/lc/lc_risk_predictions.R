# Author: Madison Coots
# Date: November 2, 2023
# ----------------------
# This file makes lung cancer risk predictions.
# Coots et al. (2024)

library(tidyverse)
library(readr)
library(lcrisks)
library(nnet) # for multinomial regression

read_path <- "~/Documents/harvard/research/lung_cancer/data/processed/lc_data.rds"
data_object_write_path <- "~/Documents/harvard/research/lung_cancer/data/processed/"
weights_path <- here::here(dirname(rstudioapi::getActiveDocumentContext()$path)) %>% 
  str_remove("code/lc") %>%
  str_c("data/lc/processed/lc_data_weights.rds")

data <- readRDS(read_path)
weights <- readRDS(weights_path)

lcrisks_data <- data %>%
  select(age,
         female,
         years_smoked,
         years_quit,
         n_cig_per_day,
         race_code,
         lung_disease,
         num_relatives_with_lc,
         bmi,
         highest_educ_level,
         personal_cancer_history,
         hypertension,
         chd,
         angina,
         heart_attack,
         other_heart_disease,
         stroke,
         diab,
         bron,
         kidney,
         liver,
         spec_eq,
         year)

# ======================================================================
# Predicting race-aware lung cancer risk (standard)
# ======================================================================

race_aware_risk <- lcrisk(lcrisks_data, nyears=5)
race_aware_lcrat <- race_aware_risk$`Number with lung cancer diagnosed per 1000 (LCRAT)` / 1000 

all_data_with_race_aware_risk <- data %>%
  mutate(race_aware_risk = race_aware_lcrat)

# ======================================================================
# Predicting race-blind 5-year lung cancer risk (non-standard)
# ======================================================================

# ====== Need to train multinomial model to compute P(R = r | X) =======

# Omitted angina, other_heart_disease, kidney, liver, spec_eq 
# because we do not have this info in NLST
# Omitted year because it is the same for all obs. in NLST
race_weight_model_formula <- race_str ~ age + female + years_smoked +
  years_quit + n_cig_per_day + lung_disease + num_relatives_with_lc + bmi +
  highest_educ_level + personal_cancer_history + hypertension + chd + 
  heart_attack + stroke + diab + bron

race_weight_model <- multinom(race_weight_model_formula,
                              data = data)

race_weights <- predict(race_weight_model, data, type = "probs") %>%
  data.frame() %>%
  rename(white_weight = White,
         black_weight = Black,
         hispanic_weight = Hispanic,
         asian_weight = Asian.PI) %>%
  mutate(pid = data$pid) %>% 
  drop_na()

# ======================================================================

white_risk <- lcrisk(lcrisks_data %>% mutate(race_code = 0), nyears = 5)
white_risk <- white_risk$`Number with lung cancer diagnosed per 1000 (LCRAT)` / 1000 
black_risk <- lcrisk(lcrisks_data %>% mutate(race_code = 1), nyears = 5)
black_risk <- black_risk$`Number with lung cancer diagnosed per 1000 (LCRAT)` / 1000 
hispanic_risk <- lcrisk(lcrisks_data %>% mutate(race_code = 2), nyears = 5)
hispanic_risk <- hispanic_risk$`Number with lung cancer diagnosed per 1000 (LCRAT)` / 1000 
asian_risk <- lcrisk(lcrisks_data %>% mutate(race_code = 3), nyears = 5)
asian_risk <- asian_risk$`Number with lung cancer diagnosed per 1000 (LCRAT)` / 1000 

all_data_with_race_blind_risk <- data %>%
  mutate(white_risk = white_risk,
         black_risk = black_risk,
         hispanic_risk = hispanic_risk,
         asian_risk = asian_risk) %>%
  inner_join(race_weights, by = c("pid")) %>% # inner join avoids NA's
  mutate(race_blind_risk = (
    white_weight * white_risk +
      black_weight * black_risk +
      hispanic_weight * hispanic_risk +
      asian_weight * asian_risk
  )) %>%
  select(-white_risk, -black_risk, -hispanic_risk, -asian_risk,
         -white_weight, -black_weight, -hispanic_weight, -asian_weight)

# ======================================================================
# Saving data and both predictions (and attaching weights)
# ======================================================================

all_data_with_both_risk <-
  all_data_with_race_aware_risk %>%
  inner_join(all_data_with_race_blind_risk %>% # inner join avoids NA's
              select(pid, race_blind_risk), by = c("pid")) %>%
  inner_join(weights, by = c("pid"))

saveRDS(all_data_with_both_risk, file = paste(data_object_write_path, "all_data_with_lc_risk.rds", sep = ""))

