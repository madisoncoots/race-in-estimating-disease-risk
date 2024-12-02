# Author: Madison Coots
# Date: September 27, 2023
# ----------------------
# This file makes cardiovascular risk predictions.
# Coots et al. (2023)

library(tidyverse)
library(readr)
library(stringr)
library(here)
library(rstudioapi)

working_directory <- getwd()
if (basename(working_directory) == "code") {
  # being called from master script
  directory_path <- here(working_directory, 'cvd/')
} else {
  # being run on its own
  directory_path <- dirname(getActiveDocumentContext()$path)
}

source(here(directory_path, "2013_pce_constants.R")) # For PCE group-specific coefficients

data_object_read_path <- directory_path %>% str_remove("code/cvd") %>% str_c("data/cvd/processed/cvd_data.rds")
data_object_write_path <- directory_path %>% str_remove("code/cvd") %>% str_c("data/cvd/processed/")

data <- readRDS(data_object_read_path)

# ======================================================================
# Function that constructs the features used to compute 10-yr ASCVD risk
# ======================================================================
make_model_data <- function(data) {
  model_data <- data %>%
    mutate(ln_age = log(ridageyr),
           ln_age_sq = ln_age^2,
           ln_tchol = log(lbxtc),
           ln_age_ln_tchol = ln_age * ln_tchol,
           ln_hdl = log(lbdhdd),
           ln_age_ln_hdl = ln_age * ln_hdl,
           ln_treated_sys_bp = if_else(
             hypertension_treatment, log(sys_bp), 0
           ),
           ln_age_ln_treated_sys_bp = ln_age * ln_treated_sys_bp,
           ln_untreated_sys_bp = if_else(
             !hypertension_treatment, log(sys_bp), 0
           ),
           ln_age_ln_untreated_sys_bp = ln_age * ln_untreated_sys_bp,
           ln_age_smokes = ln_age * smokes)  %>%
    mutate_if(is.logical, as.integer)
  return(model_data)
}

model_data <- make_model_data(data)

# ======================================================================
# Function that predicts 10-yr ASCVD risk given input data
# ======================================================================
predict_10_yr_cvd_risk <- function(coef, data, mean_val, baseline_survival) {
  covariates <- names(coef)
  data <- data %>% 
    select(any_of(covariates))
  formula <- reformulate(paste(c("-1 + ", paste(names(coef), collapse = " + ")), collapse = " "))
  model_data <- model.matrix(formula, data)[,]
  dot_product <- c(model_data %*% coef)
  ten_yr_risk <- 1 - (baseline_survival^exp(dot_product - mean_val))
  return(ten_yr_risk)
}

# ======================================================================
# Test cases
# ======================================================================

test_observation <- make_model_data(
  data.frame(ridageyr = 55,
             lbxtc = 213,
             lbdhdd = 50,
             sys_bp = 120,
             hypertension_treatment = FALSE,
             smokes= FALSE,
             diabetes = FALSE)
)

# White women (should be ~2.1%)
predict_10_yr_cvd_risk(coef = white_women_coef,
                       data = test_observation,
                       mean_val = white_women_mean_val,
                       baseline_survival = white_women_baseline_survival)

# Black women (should be ~3.0%)
predict_10_yr_cvd_risk(coef = black_women_coef,
                       data = test_observation,
                       mean_val = black_women_mean_val,
                       baseline_survival = black_women_baseline_survival)

# White men (should be ~5.3%)
predict_10_yr_cvd_risk(coef = white_men_coef,
                       data = test_observation,
                       mean_val = white_men_mean_val,
                       baseline_survival = white_men_baseline_survival)

# Black men (should be ~6.1%)
predict_10_yr_cvd_risk(coef = black_men_coef,
                       data = test_observation,
                       mean_val = black_men_mean_val,
                       baseline_survival = black_men_baseline_survival)


# ======================================================================
# Making race-aware ASCVD predictions (standard)
# ======================================================================

white_women_data <- model_data %>%
  filter(race == "White", gender == "Woman") %>%
  mutate(race_aware_ascvd_risk = predict_10_yr_cvd_risk(coef = white_women_coef,
                                                        data = .,
                                                        mean_val = white_women_mean_val,
                                                        baseline_survival = white_women_baseline_survival))

black_women_data <- model_data %>%
  filter(race == "Black", gender == "Woman")  %>%
  mutate(race_aware_ascvd_risk = predict_10_yr_cvd_risk(coef = black_women_coef,
                                                        data = .,
                                                        mean_val = black_women_mean_val,
                                                        baseline_survival = black_women_baseline_survival))

asian_women_data <- model_data %>%
  filter(race == "Asian", gender == "Woman") %>%
  mutate(race_aware_ascvd_risk = predict_10_yr_cvd_risk(coef = white_women_coef,
                                                        data = .,
                                                        mean_val = white_women_mean_val,
                                                        baseline_survival = white_women_baseline_survival))

hispanic_women_data <- model_data %>%
  filter(race == "Hispanic", gender == "Woman") %>%
  mutate(race_aware_ascvd_risk = predict_10_yr_cvd_risk(coef = white_women_coef,
                                                        data = .,
                                                        mean_val = white_women_mean_val,
                                                        baseline_survival = white_women_baseline_survival))

white_men_data <- model_data %>%
  filter(race == "White", gender == "Man")  %>%
  mutate(race_aware_ascvd_risk = predict_10_yr_cvd_risk(coef = white_men_coef,
                                                        data = .,
                                                        mean_val = white_men_mean_val,
                                                        baseline_survival = white_men_baseline_survival))

black_men_data <- model_data %>%
  filter(race == "Black", gender == "Man")  %>%
  mutate(race_aware_ascvd_risk = predict_10_yr_cvd_risk(coef = black_men_coef,
                                                        data = .,
                                                        mean_val = black_men_mean_val,
                                                        baseline_survival = black_men_baseline_survival))

asian_men_data <- model_data %>%
  filter(race == "Asian", gender == "Man")  %>%
  mutate(race_aware_ascvd_risk = predict_10_yr_cvd_risk(coef = white_men_coef,
                                                        data = .,
                                                        mean_val = white_men_mean_val,
                                                        baseline_survival = white_men_baseline_survival))

hispanic_men_data <- model_data %>%
  filter(race == "Hispanic", gender == "Man")  %>%
  mutate(race_aware_ascvd_risk = predict_10_yr_cvd_risk(coef = white_men_coef,
                                                        data = .,
                                                        mean_val = white_men_mean_val,
                                                        baseline_survival = white_men_baseline_survival))

all_data_with_race_aware_ascvd <- bind_rows(
  white_women_data,
  black_women_data,
  asian_women_data,
  hispanic_women_data,
  white_men_data,
  black_men_data,
  asian_men_data,
  hispanic_men_data,
)


# ======================================================================
# Making race-blind ASCVD predictions (non-standard)
# ======================================================================

# ====== Need to train multinomial model to compute P(R = r | X) =======
race_weight_model_formula <- is_black ~ gender + ln_age + ln_age_sq + ln_tchol + 
  ln_age_ln_tchol + ln_hdl + ln_age_ln_hdl + ln_treated_sys_bp + 
  ln_age_ln_treated_sys_bp + ln_untreated_sys_bp + ln_age_ln_untreated_sys_bp + 
  ln_age_smokes

race_weight_model_data <- model_data %>% 
  mutate(is_black = race == "Black")

race_weight_model <- glm(race_weight_model_formula,
    data = race_weight_model_data,
    family = "binomial",
    weights = round(wtmec8yr/1000)) # glm complains when weights aren't ints

race_weights <- predict(race_weight_model, model_data, type = "response") %>%
  data.frame() %>%
  rename(black_weight = ".") %>%
  mutate(not_black_weight = 1 - black_weight) %>%
  mutate(seqn = model_data$seqn)
# ======================================================================

women_data <- model_data %>%
  filter(gender == "Woman") %>%
  mutate(not_black_ascvd_risk = 
           predict_10_yr_cvd_risk(coef = white_women_coef,
                                  data = .,
                                  mean_val = white_women_mean_val,
                                  baseline_survival = white_women_baseline_survival),
         black_ascvd_risk = 
           predict_10_yr_cvd_risk(coef = black_women_coef,
                                  data = .,
                                  mean_val = black_women_mean_val,
                                  baseline_survival = black_women_baseline_survival)) %>%
  left_join(race_weights, by = c("seqn")) %>%
  mutate(race_blind_ascvd_risk = (not_black_weight * not_black_ascvd_risk + 
                                    black_weight * black_ascvd_risk)) %>%
  select(-not_black_ascvd_risk, -black_ascvd_risk, -not_black_weight, -black_weight)

men_data <- model_data %>%
  filter(gender == "Man") %>%
  mutate(not_black_ascvd_risk = 
           predict_10_yr_cvd_risk(coef = white_men_coef,
                                  data = .,
                                  mean_val = white_men_mean_val,
                                  baseline_survival = white_men_baseline_survival),
         black_ascvd_risk = 
           predict_10_yr_cvd_risk(coef = black_men_coef,
                                  data = .,
                                  mean_val = black_men_mean_val,
                                  baseline_survival = black_men_baseline_survival)) %>%
  left_join(race_weights, by = c("seqn")) %>% 
  mutate(race_blind_ascvd_risk = (not_black_weight * not_black_ascvd_risk + 
                                    black_weight * black_ascvd_risk)) %>%
  select(-not_black_ascvd_risk, -black_ascvd_risk, -not_black_weight, -black_weight)

all_data_with_race_blind_ascvd <- bind_rows(
  women_data,
  men_data
)  

# ======================================================================
# Saving data and both predictions
# ======================================================================

all_data_with_both_ascvd <-
  all_data_with_race_aware_ascvd %>%
  left_join(all_data_with_race_blind_ascvd %>%
              select(seqn, race_blind_ascvd_risk),
            by = c("seqn"))
  
saveRDS(all_data_with_both_ascvd, file = paste(data_object_write_path, "cvd_data_with_risk.rds", sep = ""))  
  