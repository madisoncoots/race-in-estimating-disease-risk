# Author: Madison Coots
# Date: November 30, 2024
# ======================
# Code for all breast cancer figures, tables, and statistics.
# Coots et al. (2024)

library(tidyverse)
library(readr)
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

source(here(directory_path, "/colors.R")) # For figure color maps
data_object_read_path <- directory_path %>% str_remove("code/bc") %>% str_c("data/bc/processed/bc_data_with_risk.rds")
data <- readRDS(data_object_read_path) %>%
  rename(race = race_str)

theme_set(theme_bw(base_size = 15))

################################## Parameters ##################################

risk_score_upper_bound <- 0.025
incidence_upper_bound <- 0.026
screening_thresh <- 0.0166
optimal_thresh <- screening_thresh
screening_cost <- -screening_thresh
utility_reward <- 1
per_capita_n <- 10000

################################### Figures ####################################

################################################################################
# Figure 2a: Race-blind calibration plot

bc_race_blind_calib <- data %>%
    select(race, race_blind_risk, race_aware_risk, wtmec8yr) %>%
    mutate(bin = floor((race_blind_risk  + 0.0025) * 100 * 2) / 2 / 100,
           bucket = ntile(race_blind_risk / 0.025, 10)) %>%
    group_by(race, bin) %>%
    summarize(n_in_bin = sum(wtmec8yr),
              bin_avg = sum(race_blind_risk * wtmec8yr) / sum(wtmec8yr),
              cvd_prev = sum(race_aware_risk * wtmec8yr) / sum(wtmec8yr)) %>%
    mutate(disease = "Breast cancer") %>% # For eventual faceting
    ggplot(aes(x=bin_avg, y=cvd_prev, color=race)) +
    geom_vline(xintercept=screening_thresh) +
    geom_line() +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    xlab("Race-blind BC risk") +
    ylab("True BC Risk") +
    scale_y_continuous(labels = scales::percent,
                       breaks = seq(0.0, incidence_upper_bound, 0.005)) +
    scale_x_continuous(labels = scales::percent,
                       breaks = seq(0.0, risk_score_upper_bound, 0.005)) +
    coord_cartesian(xlim = c(0, risk_score_upper_bound), ylim = c(0, incidence_upper_bound)) +
    theme(legend.title = element_blank(),
          legend.position = c(0.18, 0.81)) +
    scale_color_manual(values=group_color_map,
                       breaks =group_names)

################################################################################
# Figure 2b: Race-blind scatter plot

bc_race_blind_scatter <- data %>%
  select(race, race_blind_risk, race_aware_risk, wtmec8yr) %>%
  ggplot(aes(x=race_blind_risk, y=race_aware_risk, color=race)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(shape = 1) +
  annotate("rect", xmin = -1, xmax = screening_thresh, ymin = -1, ymax = screening_thresh,
           alpha = 0.6, fill="lightgray") +
  annotate("rect", xmin = screening_thresh, xmax = 1, ymin = screening_thresh, ymax = 1,
           alpha = 0.6, fill="lightgray") +
  geom_vline(xintercept=screening_thresh) +
  geom_hline(yintercept=screening_thresh) +
  xlab("Race-blind BC risk") +
  ylab("True BC risk") +
  scale_color_manual(values=group_color_map,
                     breaks =group_names) +
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.86)) +
  coord_cartesian(xlim = c(0, risk_score_upper_bound), ylim = c(0, incidence_upper_bound)) +
  scale_x_continuous(labels = scales::percent,
                     breaks = seq(0.0, risk_score_upper_bound, 0.005)) +
  scale_y_continuous(labels = scales::percent,
                     breaks = seq(0.0, incidence_upper_bound, 0.005))

################################################################################
# Figure 3a: Utility dot plot

utility_gains_by_group <- 
  data %>% 
  mutate(p_disease = race_aware_risk,
         race_aware_risk_score = race_aware_risk,
         race_blind_risk_score = race_blind_risk,
         race_aware_decision = race_aware_risk_score > optimal_thresh,
         race_blind_decision = race_blind_risk_score > optimal_thresh,
         expected_race_aware_utility = (screening_cost + utility_reward * p_disease) * race_aware_decision,
         expected_race_blind_utility = (screening_cost + utility_reward * p_disease) * race_blind_decision,
         utility_difference = expected_race_aware_utility - expected_race_blind_utility
  ) %>%
  drop_na(utility_difference) %>% # drop observations that have an NA prediction due to missing data
  group_by(race) %>%
  summarize(average_utility_diff = sum(utility_difference * wtmec8yr)/sum(wtmec8yr) * per_capita_n,
            n = sum(wtmec8yr)) %>%
  arrange(desc(average_utility_diff))

utility_gains_overall <- utility_gains_by_group %>%
  group_by() %>%
  summarize(average_utility_diff = sum(average_utility_diff * n) / sum(n)) %>%
  pull(average_utility_diff)

dot_order <- utility_gains_by_group %>% pull(race)

bc_utility_dot_plot <- utility_gains_by_group %>%
  mutate(disease = "Breast cancer") %>% # For easier faceting
  ggplot(aes(x=average_utility_diff, y=factor(race, levels=dot_order), color=race, fill = race)) +
  geom_point(size=4) +
  geom_bar(stat="identity", width = 0.03) +
  geom_vline(xintercept=utility_gains_overall) + 
  theme(legend.position="none") +
  xlab("Additional cases detected (at no cost) with race-aware predictions per 10,000 individuals") +
  ylab("") +
  scale_color_manual(values=group_color_map,
                     breaks =group_names) +
  scale_fill_manual(values=group_color_map,
                    breaks =group_names) +
  scale_y_discrete(labels = NULL,
                   breaks = NULL)

################################################################################
# Figure 3b: Pct. changed decisions dot plot

pct_changed_by_group <- data %>%
  select(ID, race, race_blind_risk, race_aware_risk, wtmec8yr) %>%
  mutate(race_blind_screening_decision = race_blind_risk > screening_thresh,
         race_aware_screening_decision = race_aware_risk > screening_thresh,
         changed_decision = (race_blind_screening_decision != race_aware_screening_decision)) %>%
  group_by(race) %>%
  summarize(pct_changed_decision = sum(changed_decision * wtmec8yr) / sum(wtmec8yr),
            n = sum(wtmec8yr)) %>%
  arrange(pct_changed_decision)

pct_changed_overall <- pct_changed_by_group %>%
  group_by() %>%
  summarize(pct_changed_decision = sum(pct_changed_decision * n) / sum(n)) %>%
  pull(pct_changed_decision)

dot_order <- pct_changed_by_group %>% pull(race)

bc_dot_plot <- pct_changed_by_group %>%
  mutate(disease = "Breast cancer") %>% # For easier faceting
  ggplot(aes(x=pct_changed_decision, y=factor(race, levels=dot_order), color=race, fill = race)) +
  geom_point(size=4) +
  geom_bar(stat="identity", width = 0.03) +
  geom_vline(xintercept=pct_changed_overall) + 
  theme(legend.position="none") +
  xlab("Pct. of individuals with\nchanged decisions") +
  ylab("") +
  scale_color_manual(values=group_color_map,
                     breaks =group_names) +
  scale_fill_manual(values=group_color_map,
                    breaks =group_names) +
  scale_x_continuous(labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(xlim = c(0, 1))

################################################################################
# Figure 4a: Utility under scarcity plot

# Computes scarcity threshold when we limit to the top K%
compute_blind_scarcity_thresh_using_pct <- function(data, pct_ppl_we_can_screen) {
  data %>%
    mutate(risk_score = race_blind_risk) %>%
    drop_na(risk_score) %>%
    arrange(desc(risk_score)) %>%
    mutate(cumulative_sum = cumsum(wtmec8yr),
           cumulative_pct = cumulative_sum / sum(wtmec8yr) * 100) %>%
    filter(cumulative_pct >= pct_ppl_we_can_screen) %>%
    pull(risk_score) %>%
    first(default = NA)
}

compute_aware_scarcity_thresh_using_pct <- function(data, pct_ppl_we_can_screen) {
  data %>%
    mutate(risk_score = race_aware_risk) %>%
    drop_na(risk_score) %>%
    arrange(desc(risk_score)) %>%
    mutate(cumulative_sum = cumsum(wtmec8yr),
           cumulative_pct = cumulative_sum / sum(wtmec8yr) * 100) %>%
    filter(cumulative_pct >= pct_ppl_we_can_screen) %>%
    pull(risk_score) %>%
    first(default = NA)
}

capacities <- seq(0, 100, 1)
scarcity_utility_plot_data <- c()
for (c in capacities) {
  race_blind_scarcity_thresh <- compute_blind_scarcity_thresh_using_pct(data, c)
  race_aware_scarcity_thresh <- compute_aware_scarcity_thresh_using_pct(data, c)
  avg_utility_diff_by_group <- data %>%
    mutate(p_disease = race_aware_risk,
           race_aware_risk_score = race_aware_risk,
           race_blind_risk_score = race_blind_risk,
           race_aware_decision = race_aware_risk_score > max(screening_thresh, race_aware_scarcity_thresh),
           race_blind_decision = race_blind_risk_score > max(screening_thresh, race_blind_scarcity_thresh),
           expected_race_aware_utility = (screening_cost + utility_reward * p_disease) * race_aware_decision,
           expected_race_blind_utility = (screening_cost + utility_reward * p_disease) * race_blind_decision,
           utility_difference = expected_race_aware_utility - expected_race_blind_utility
    ) %>%
    drop_na(utility_difference) %>%
    group_by(race) %>%
    summarize(average_utility_diff_dollars = sum(utility_difference * wtmec8yr) / sum(wtmec8yr) * per_capita_n) %>%
    arrange(desc(average_utility_diff_dollars)) %>%
    mutate(capacity = c,
           race_blind_scarcity_thresh = race_blind_scarcity_thresh,
           race_aware_scarcity_thresh = race_aware_scarcity_thresh)
  scarcity_utility_plot_data <- rbind(scarcity_utility_plot_data, avg_utility_diff_by_group)
}

bc_scarcity_utility_plot <- scarcity_utility_plot_data %>%
  filter(race_blind_scarcity_thresh > screening_thresh) %>%
  mutate(disease = "Breast cancer") %>%
  ggplot(aes(x = race_blind_scarcity_thresh, y = average_utility_diff_dollars, color = race)) +
  geom_hline(yintercept=0) +
  # geom_point() +
  geom_line() +
  geom_vline(xintercept=screening_thresh) +
  ylab("Gains in net benefit from the\nswitch to race-aware predictions\nper 10,000 individuals") +
  xlab("Screening threshold") +
  scale_x_continuous(labels = scales::percent) +
  theme(legend.title = element_blank(),
        legend.position =  c(0.85, 0.87)) +
  scale_color_manual(values=group_color_map,
                     breaks =group_names)

################################################################################
# Figure 4b: Pct. changed decisions under scarcity plot

capacities <- seq(0, 100, 1)
scarcity_pct_changed_plot_data <- c()
for (c in capacities) {
  race_blind_scarcity_thresh <- compute_blind_scarcity_thresh_using_pct(data, c)
  race_aware_scarcity_thresh <- compute_aware_scarcity_thresh_using_pct(data, c)
  pct_changed_by_group <- data %>%
    drop_na(race_blind_risk, race_aware_risk) %>%
    mutate(race_blind_screening_decision = race_blind_risk > max(screening_thresh, race_blind_scarcity_thresh),
           race_aware_screening_decision = race_aware_risk > max(screening_thresh, race_aware_scarcity_thresh),
           changed_decision = (race_blind_screening_decision != race_aware_screening_decision)) %>%
    group_by(race) %>%
    summarize(pct_changed_decision = sum(changed_decision * wtmec8yr) / sum(wtmec8yr),
              n = sum(wtmec8yr)) %>%
    arrange(pct_changed_decision) %>%
    mutate(capacity = c,
           race_blind_scarcity_thresh = race_blind_scarcity_thresh,
           race_aware_scarcity_thresh = race_aware_scarcity_thresh)
  scarcity_pct_changed_plot_data <- rbind(scarcity_pct_changed_plot_data, pct_changed_by_group)
}
  
bc_scarcity_pct_changed_plot <- scarcity_pct_changed_plot_data %>%
  filter(race_blind_scarcity_thresh > screening_thresh) %>%
  mutate(disease = "Breast cancer") %>%
  ggplot(aes(x = race_blind_scarcity_thresh, y = pct_changed_decision, color = race)) +
  geom_hline(yintercept=0) +
  geom_line() +
  geom_vline(xintercept=screening_thresh) +
  ylab("Pct. of individuals with changed decisions") +
  xlab("Screening threshold") +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.title = element_blank(),
        legend.position =  c(0.85, 0.87)) +
  scale_color_manual(values=group_color_map,
                     breaks =group_names) +
  coord_cartesian(ylim = c(0,1))

################################################################################
# Appendix Figure 1: Baseline utility plot

baseline_utility_by_group <- data %>% 
  mutate(p_disease = race_aware_risk,
         race_blind_risk_score = race_blind_risk,
         race_blind_decision = race_blind_risk_score > optimal_thresh,
         baseline_utility = (screening_cost + utility_reward * p_disease) * race_blind_decision,
  ) %>%
  drop_na(baseline_utility) %>% # drop observations that have an NA prediction due to missing data
  group_by(race) %>%
  summarize(average_baseline_utility = sum(baseline_utility * wtmec8yr)/sum(wtmec8yr) * per_capita_n,
            n = sum(wtmec8yr))
  
bc_baseline_utility_plot <- baseline_utility_by_group %>%
  mutate(disease = "Breast cancer") %>% # For easier faceting
  ggplot(aes(x=average_baseline_utility, y=factor(race), color=race, fill = race)) +
  geom_point(size=4) +
  geom_bar(stat="identity", width = 0.03) +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Breast cancer") +
  xlab("Cases detected with race-blind\npredictions per 10,000 individuals") +
  ylab("") +
  scale_color_manual(values=group_color_map,
                     breaks =group_names) +
  scale_fill_manual(values=group_color_map,
                    breaks =group_names) +
  scale_y_discrete(labels = NULL,
                   breaks = NULL)

baseline_utility_overall <- function(data, screening_thresh, per_capita_n) {
  bc_baseline_utility_overall <- baseline_utility(data, screening_thresh, per_capita_n) %>%
    summarise(baseline_utility_overall = sum(average_baseline_utility * n) / sum(n))
  return(bc_baseline_utility_overall)
}

################################################################################
# Appendix Figure 2: Subgroup histogram plot

subgroup_aware_data <- data %>%
  mutate(model = "Breast cancer\n\nRace-aware") %>%
  mutate(decision = race_aware_risk > screening_thresh,
         detection = decision * race_aware_risk) %>% # aware is "true" risk
  filter(detection > screening_thresh)

subgroup_blind_data <- data %>%
  mutate(model = "Breast cancer\n\nRace-unaware") %>%
  mutate(decision = race_blind_risk > screening_thresh,
         detection = decision * race_aware_risk) %>% # aware is "true" risk
  filter(detection > screening_thresh)

subgroup_data <- rbind(subgroup_aware_data,
                       subgroup_blind_data)

subgroup_means <- subgroup_data %>%
  group_by(model, race) %>%
  summarize(mean_age = sum(T1 * wtmec8yr) / sum(wtmec8yr))

bc_subgroup_histogram <- 
  subgroup_data %>%
  mutate(disease = "Breast cancer") %>%
  ggplot(aes(x = T1, y = ..density.., color = race, weight=wtmec8yr)) +
  facet_grid(model ~ race) +
  geom_histogram(binwidth = 1, fill="white", position="dodge") +
  ylab("Fraction within group") +
  xlab("Age") +
  geom_vline(data = subgroup_means, aes(xintercept =mean_age)) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = group_color_map,
                     breaks = group_names) +
  theme(legend.title = element_blank(),
        legend.position = "right")

################################################################################
# Appendix Figure 3: Sensitivity pareto plot

benefit_range <- seq(0.5, 2, 0.1)
all_results <- c()

for (k in benefit_range) {
  result <- 
    data %>%
    mutate(new_threshold = optimal_thresh / k,
           aware_screened = race_aware_risk > new_threshold,
           detection = aware_screened * race_aware_risk) %>% 
    group_by(race) %>%
    summarize(num_aware_screened = sum(aware_screened * wtmec8yr), 
              frac_aware_screened = num_aware_screened / sum(wtmec8yr),
              num_aware_detections = sum(detection * wtmec8yr),
              frac_aware_detections = num_aware_detections / sum(wtmec8yr * race_aware_risk),
              new_threshold = optimal_thresh / k,
              disease = "Breast cancer")
  all_results <- rbind(all_results, result)
}

status_quo_aware <- all_results %>%
  filter(new_threshold == optimal_thresh) %>% 
  mutate(model="Race-aware")

status_quo_blind <- data %>%
  mutate(disease = "Breast cancer",
         blind_screened = race_blind_risk > optimal_thresh,
         detection = blind_screened * race_aware_risk) %>%
  group_by(race) %>%
  summarize(num_blind_screened = sum(blind_screened * wtmec8yr), 
            frac_blind_screened = num_blind_screened / sum(wtmec8yr),
            num_blind_detections = sum(detection * wtmec8yr),
            frac_blind_detections = num_blind_detections / sum(wtmec8yr * race_aware_risk),
            new_threshold = optimal_thresh,
            disease = "Breast cancer") %>% 
  mutate(model="Race-unaware")

bc_sensitivity_pareto_plot <- 
  all_results %>%
  ggplot(aes(x=frac_aware_screened, y=frac_aware_detections, color = race)) +
  facet_grid(disease ~ race) + 
  geom_line() +
  geom_point(data = status_quo_aware, aes(x=frac_aware_screened, y = frac_aware_detections, shape = factor(model)), size = 3) +
  geom_point(data = status_quo_blind, aes(x=frac_blind_screened, y = frac_blind_detections, shape = factor(model)), size = 3) +
  ylab("Recall") +
  xlab("Fraction recommended chemoprevention") +
  scale_x_continuous(labels = scales::percent,
                     breaks = seq(0.0, 1, 0.2)) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = group_color_map,
                     breaks = group_names) +
  theme(legend.title = element_blank(),
        legend.position = "right") +
  scale_shape_manual(name = "", 
                     values = c("Race-aware" = 19, "Race-unaware" = 1)) +
  guides(shape = guide_legend(override.aes = list(fill = c('Race-aware' = "white", 'Race-unaware' = 'black'),
                                                  shape = c('Race-aware' = 19, 'Race-unaware' = 1))))

################################## Statistics ##################################

################################################################################
# Utility gains by group

utility_gains_by_group

################################################################################
# Utility gain whole population

utility_gains_overall

################################################################################
# Pct. of people with the same screening decision under both models by group

pct_same_by_group <- data %>%
  select(ID, race_blind_risk, race_aware_risk, wtmec8yr, race) %>%
  mutate(race_blind_screening_decision = race_blind_risk > screening_thresh,
         race_aware_screening_decision = race_aware_risk > screening_thresh,
         same_decision = (race_blind_screening_decision == race_aware_screening_decision)) %>%
  group_by(race) %>%
  summarize(pct_blind_decision = sum(race_blind_screening_decision * wtmec8yr) / sum(wtmec8yr) * 100,
            pct_aware_decision = sum(race_aware_screening_decision * wtmec8yr) / sum(wtmec8yr) * 100,
            pct_same_decision = sum(same_decision * wtmec8yr) / sum(wtmec8yr) * 100) %>%
  select(pct_same_decision)

pct_same_by_group

################################################################################
# Pct. of people with the same screening decision under both models whole pop

pct_same_overall <- data %>%
  select(ID, race_blind_risk, race_aware_risk, wtmec8yr) %>%
  mutate(race_blind_screening_decision = race_blind_risk > screening_thresh,
         race_aware_screening_decision = race_aware_risk > screening_thresh,
         same_decision = (race_blind_screening_decision == race_aware_screening_decision)) %>%
  summarize(pct_blind_decision = sum(race_blind_screening_decision * wtmec8yr) / sum(wtmec8yr) * 100,
            pct_aware_decision = sum(race_aware_screening_decision * wtmec8yr) / sum(wtmec8yr) * 100,
            pct_same_decision = sum(same_decision * wtmec8yr) / sum(wtmec8yr) * 100) %>%
  select(pct_same_decision)

pct_same_overall

################################################################################
# Baseline utility by group

baseline_utility_by_group

################################################################################
# Baseline utility whole pop

baseline_utility_overall <- baseline_utility_by_group %>%
  summarise(baseline_utility_overall = sum(average_baseline_utility * n) / sum(n))

baseline_utility_overall

################################################################################
# Appendix Table: Verifying optimal thresholds under both models

compute_pop_utility <- function(threshold, 
                                screening_cost, 
                                utility_reward,
                                risk_scores, 
                                p_disease, 
                                weights) {
  data <- data.frame(
    risk_score = risk_scores,
    p_disease = p_disease
  )
  population_utility <- data %>%
    mutate(screening_decision = risk_score >= threshold,
           utility = (screening_cost + (utility_reward) * p_disease) * screening_decision
    ) %>%
    summarize(average_utility = sum(utility * weights)/sum(weights)) %>%
    pull()
  return(population_utility)
}

race_aware_optimal_thresh <- optimize(compute_pop_utility, 
                                      interval = c(0,0.05),
                                      screening_cost = screening_cost,
                                      utility_reward = utility_reward,
                                      risk_scores = data$race_aware_risk,
                                      p_disease = data$race_aware_risk,
                                      weights = data$wtmec8yr,
                                      maximum = TRUE)[1]$maximum

race_blind_optimal_thresh <- optimize(compute_pop_utility, 
                                      interval = c(0,0.05),
                                      screening_cost = screening_cost,
                                      utility_reward = utility_reward,
                                      risk_scores = data$race_blind_risk,
                                      p_disease = data$race_aware_risk,
                                      weights = data$wtmec8yr,
                                      maximum = TRUE)[1]$maximum

threshold_table <- data.frame(
  disease = "Breast cancer",
  recommended_thresh = optimal_thresh,
  race_aware_optimal_thresh = race_aware_optimal_thresh,
  race_blind_optimal_thresh = race_blind_optimal_thresh
)

threshold_table
