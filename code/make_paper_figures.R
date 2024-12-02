# Author: Madison Coots
# Date: December 2, 2024
# ======================
# Code for all figures. 
# Coots et al. (2024)

library(tidyverse)
library(rstudioapi)
library(here)
library(patchwork)

options(scipen = 999) # Avoid scientific notation

directory_path <- dirname(getActiveDocumentContext()$path)
setwd(directory_path)
figure_save_path <- getwd() %>% str_remove("code") %>% str_c("figures/")

cvd_script_path <- here(directory_path, "cvd/figures_and_analysis.R")
bc_script_path <- here(directory_path, "bc/figures_and_analysis.R")
lc_figures_path <- directory_path %>% str_remove("code") %>% str_c("data/lc/processed/figure_objects")

# =============================================================================
# Run figures_and_analysis.R scripts for each disease to load the figures into
# the R studio environment.
# =============================================================================
source(cvd_script_path)
source(bc_script_path)

# =============================================================================
# Load pre-generated figures for lung cancer because NLST data is not publicly
# releasable and so figures cannot be generated directly.
# =============================================================================

lc_race_blind_calib <- readRDS(path(lc_figures_path, "lc_race_blind_calib.rds"))
lc_race_blind_scatter <- readRDS(path(lc_figures_path, "lc_race_blind_scatter.rds"))
lc_utility_dot_plot <- readRDS(path(lc_figures_path, "lc_utility_dot_plot.rds"))
lc_dot_plot <- readRDS(path(lc_figures_path, "lc_dot_plot.rds"))
lc_scarcity_utility_plot <- readRDS(path(lc_figures_path, "lc_scarcity_utility_plot.rds"))
lc_scarcity_pct_changed_plot <- readRDS(path(lc_figures_path, "lc_scarcity_pct_changed_plot.rds"))
lc_baseline_utility_plot <- readRDS(path(lc_figures_path, "lc_baseline_utility_plot.rds"))
lc_subgroup_histogram <- readRDS(path(lc_figures_path, "lc_subgroup_histogram.rds"))
lc_sensitivity_pareto_plot <- readRDS(path(lc_figures_path, "lc_sensitivity_pareto_plot.rds"))

# =============================================================================
# Figure 2: Race-blind calibration and scatter plots.
# =============================================================================

# Removing legends, axis labels, and adding facet labels for better presentation 
# in a faceted grid
cvd_race_blind_calib <- cvd_race_blind_calib + 
  theme(legend.position="none")  + 
  facet_wrap(vars(fct_rev(disease)), ncol = 1) +
  xlab("") +
  ylab("Race-aware risk")

cvd_race_blind_scatter <- cvd_race_blind_scatter + 
  theme(legend.position="none") +
  xlab("Race-unaware risk") +
  ylab("Race-aware risk")

bc_race_blind_calib <- bc_race_blind_calib + 
  theme(legend.position="none") + 
  facet_wrap(vars(fct_rev(disease)), ncol = 1) +
  xlab("") +
  ylab("")

bc_race_blind_scatter <- bc_race_blind_scatter + 
  theme(legend.position="none") +
  xlab("Race-unaware risk") +
  ylab("")

lc_race_blind_calib <- lc_race_blind_calib + 
  theme(legend.position="right") + 
  xlab("") +
  ylab("") +
  facet_wrap(vars(fct_rev(disease)), ncol = 1)

lc_race_blind_scatter <- lc_race_blind_scatter + 
  guides(color = "none") + # drop the legend for the points
  theme(legend.position="right") + 
  xlab("Race-unaware risk") +
  ylab("")

calibration_scatter_grid <- cvd_race_blind_calib + 
  bc_race_blind_calib + lc_race_blind_calib +
  cvd_race_blind_scatter + 
  bc_race_blind_scatter + lc_race_blind_scatter +
  plot_layout(ncol = 3, nrow = 2)

calibration_scatter_grid

ggsave(paste(figure_save_path, "figure_2.pdf", sep = ""),
       width = 13,
       height = 7)

# =============================================================================
# Figure 3: Utility and pct. changed decisions dot plot
# =============================================================================

cvd_dot_plot_data <- ggplot_build(cvd_dot_plot)$plot$data
bc_dot_plot_data <- ggplot_build(bc_dot_plot)$plot$data
lc_dot_plot_data <- ggplot_build(lc_dot_plot)$plot$data

cvd_utility_dot_plot_data <- ggplot_build(cvd_utility_dot_plot)$plot$data
bc_utility_dot_plot_data <- ggplot_build(bc_utility_dot_plot)$plot$data
lc_utility_dot_plot_data <- ggplot_build(lc_utility_dot_plot)$plot$data

dot_plot_data_all <- rbind(cvd_dot_plot_data,
                           bc_dot_plot_data, 
                           lc_dot_plot_data) %>%
  group_by(disease) %>%
  arrange(desc(pct_changed_decision)) %>%
  mutate(pop_average = sum(pct_changed_decision * n) / sum(n),
         order = row_number(),
         disease = factor(disease, levels=c("Cardiovascular disease",
                                            "Breast cancer", "Lung cancer"))) %>%
  arrange(disease, order)

utility_dot_plot_data_all <- rbind(cvd_utility_dot_plot_data,
                                   bc_utility_dot_plot_data, 
                                   lc_utility_dot_plot_data) %>%
  group_by(disease) %>%
  arrange(desc(average_utility_diff)) %>%
  mutate(pop_average = sum(average_utility_diff * n) / sum(n),
         order = row_number(),
         disease = factor(disease, levels=c("Cardiovascular disease",
                                            "Breast cancer", "Lung cancer"))) %>%
  arrange(disease, order)

dot_plot <- dot_plot_data_all %>% 
  ggplot(aes(x=pct_changed_decision, y=factor(-order), color=race, fill = race)) +
  facet_wrap(~disease, ncol = 4) +
  geom_point(size=4) +
  geom_bar(stat="identity", width = 0.03, show.legend = FALSE) +
  geom_vline(data = dot_plot_data_all, aes(xintercept = pop_average), 
             show.legend = FALSE) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  xlab("Fraction of individuals with changed decisions") +
  ylab("") +
  scale_color_manual(values=group_color_map,
                     breaks =group_names) +
  scale_fill_manual(values=group_color_map,
                    breaks =group_names) +
  scale_x_continuous(labels = scales::percent,
                     breaks = seq(0, 1, 0.25),
                     minor_breaks = c()) +
  scale_y_discrete(labels = NULL,
                   breaks = NULL) +
  coord_cartesian(xlim = c(0, 1))

utility_dot_plot <- utility_dot_plot_data_all %>% 
  ggplot(aes(x=average_utility_diff, y=factor(-order), color=race, fill = race)) +
  facet_wrap(~disease, ncol = 4) +
  geom_point(size=4) +
  geom_bar(stat="identity", width = 0.03, show.legend = FALSE) +
  geom_vline(data = utility_dot_plot_data_all, aes(xintercept = pop_average), 
             show.legend = FALSE) +
  theme(legend.title = element_blank(),
        legend.position = "right") + 
  xlab("Gains in net benefit from the switch to race-aware predictions per 10,000 individuals") +
  ylab("") +
  scale_color_manual(values=group_color_map,
                     breaks =group_names) +
  scale_fill_manual(values=group_color_map,
                    breaks =group_names) +
  scale_y_discrete(labels = NULL,
                   breaks = NULL) +
  coord_cartesian(xlim=c(0,30.1))

dot_plot_grid <- utility_dot_plot + dot_plot + 
  plot_layout(ncol = 1, nrow = 2)

dot_plot_grid

ggsave(paste(figure_save_path, "figure_3.pdf", sep = ""),
       width = 12,
       height = 5)

# =============================================================================
# Figure 4: Scarcity plots
# =============================================================================

cvd_scarcity_utility_plot_data <- ggplot_build(cvd_scarcity_utility_plot)$plot$data
bc_scarcity_utility_plot_data <- ggplot_build(bc_scarcity_utility_plot)$plot$data
lc_scarcity_utility_plot_data <- ggplot_build(lc_scarcity_utility_plot)$plot$data

cvd_scarcity_pct_changed_plot_data <- ggplot_build(cvd_scarcity_pct_changed_plot)$plot$data
bc_scarcity_pct_changed_plot_data <- ggplot_build(bc_scarcity_pct_changed_plot)$plot$data
lc_scarcity_pct_changed_plot_data <- ggplot_build(lc_scarcity_pct_changed_plot)$plot$data

scarcity_utility_plot_data_all <- rbind(cvd_scarcity_utility_plot_data,
                                        bc_scarcity_utility_plot_data, 
                                        lc_scarcity_utility_plot_data) %>%
  group_by(disease) %>%
  mutate(disease = factor(disease, levels=c("Cardiovascular disease",
                                            "Breast cancer", "Lung cancer")))

scarcity_pct_changed_plot_data_all <- rbind(cvd_scarcity_pct_changed_plot_data,
                                            bc_scarcity_pct_changed_plot_data,
                                            lc_scarcity_pct_changed_plot_data) %>%
  group_by(disease) %>%
  mutate(disease = factor(disease, levels=c("Cardiovascular disease",
                                            "Breast cancer", "Lung cancer")))

cvd_top <- cvd_scarcity_utility_plot + 
  facet_wrap(~disease, ncol = 3) + 
  theme(legend.title = element_blank(),
        legend.position =  "none") +
  xlab("") +
  coord_cartesian(xlim = c(0,0.8),
                  ylim = c(-25,100))

bc_top <- bc_scarcity_utility_plot + 
  facet_wrap(~disease, ncol = 3) + 
  theme(legend.title = element_blank(),
        legend.position =  "none",
        axis.text.y = element_blank(),  # Remove y-axis labels
        axis.ticks.y = element_blank()) + 
  ylab("") +
  xlab("") +
  coord_cartesian(xlim = c(0,0.04),
                  ylim = c(-25,100))

lc_top <- lc_scarcity_utility_plot + 
  facet_wrap(~disease, ncol = 3) + 
  theme(legend.title = element_blank(),
        legend.position =  "right",
        axis.text.y = element_blank(),  # Remove y-axis labels
        axis.ticks.y = element_blank()) + 
  ylab("") +
  xlab("") +
  coord_cartesian(xlim = c(0,0.3), ylim = c(-25,100))

cvd_bottom <- cvd_scarcity_pct_changed_plot + 
  facet_wrap(~disease, ncol = 3) + 
  theme(legend.title = element_blank(),
        legend.position =  "none") +
  xlab("") +
  ylab("Fraction of individuals\nwith changed decisions") +
  coord_cartesian(xlim = c(0,0.8),
                  ylim = c(0,1))

bc_bottom <- bc_scarcity_pct_changed_plot + 
  facet_wrap(~disease, ncol = 3) + 
  theme(legend.title = element_blank(),
        legend.position =  "none",
        axis.text.y = element_blank(),  # Remove y-axis labels
        axis.ticks.y = element_blank()) + 
  ylab("") +
  coord_cartesian(xlim = c(0,0.04),
                  ylim = c(0,1))

lc_bottom <- lc_scarcity_pct_changed_plot + 
  facet_wrap(~disease, ncol = 3) + 
  theme(legend.title = element_blank(),
        legend.position =  "none",
        axis.text.y = element_blank(),  # Remove y-axis labels
        axis.ticks.y = element_blank()) + 
  ylab("") +
  xlab("") +
  coord_cartesian(xlim = c(0,0.3), ylim = c(0,1))

scarcity_grid <- cvd_top + bc_top + lc_top + 
  cvd_bottom + bc_bottom + lc_bottom + 
  plot_layout(ncol = 3, nrow = 2)

scarcity_grid

ggsave(paste(figure_save_path, "figure_4.pdf", sep = ""),
       width = 13,
       height = 7)

# =============================================================================
# Appendix Figure 1: Baseline utility plot
# =============================================================================

cvd_baseline_utility_plot_data <- ggplot_build(cvd_baseline_utility_plot)$plot$data
bc_baseline_utility_plot_data <- ggplot_build(bc_baseline_utility_plot)$plot$data
lc_baseline_utility_plot_data <- ggplot_build(lc_baseline_utility_plot)$plot$data

baseline_utility_plot_data_all <- rbind(cvd_baseline_utility_plot_data,
                                        bc_baseline_utility_plot_data,
                                        lc_baseline_utility_plot_data 
)%>%
  group_by(disease) %>%
  arrange(desc(average_baseline_utility)) %>%
  mutate(pop_average = sum(average_baseline_utility * n) / sum(n),
         order = row_number(),
         disease = factor(disease, levels=c("Cardiovascular disease",
                                            "Breast cancer", "Lung cancer"))) %>%
  arrange(disease, order)

baseline_utility_plot <- baseline_utility_plot_data_all %>% 
  ggplot(aes(x=average_baseline_utility, y=factor(-order), color=race, fill = race)) +
  facet_wrap(~disease, ncol = 4) +
  geom_point(size=4) +
  geom_bar(stat="identity", width = 0.03, show.legend = FALSE) +
  theme(legend.title = element_blank(),
        legend.position = "right") + 
  xlab("Net benefit from race-blind predictions per 10,000 individuals") +
  ylab("") +
  scale_color_manual(values=group_color_map,
                     breaks =group_names) +
  scale_fill_manual(values=group_color_map,
                    breaks =group_names) +
  scale_y_discrete(labels = NULL,
                   breaks = NULL)

baseline_utility_plot

ggsave(paste(figure_save_path, "appendix_figure_1.pdf", sep = ""),
       width = 12,
       height = 3)

# =============================================================================
# Appendix Figure 2: Subgroup histogram plot
# =============================================================================

new_cvd_subgroup_histogram <- 
  cvd_subgroup_histogram +
  theme(legend.title = element_blank(),
        legend.position = "none")

new_lc_subgroup_histogram <- 
  lc_subgroup_histogram + 
  theme(legend.title = element_blank(),
        legend.position = "none")

subgroup_histogram_grid <- new_cvd_subgroup_histogram +
  bc_subgroup_histogram +
  new_lc_subgroup_histogram +
  plot_layout(ncol=1, nrow=3)

subgroup_histogram_grid

ggsave(paste(figure_save_path, "appendix_figure_2.pdf", sep = ""),
       width = 12,
       height = 15)

# =============================================================================
# Appendix Figure 3: Sensitivity pareto plot
# =============================================================================

new_cvd_sensitivity_pareto_plot <-
  cvd_sensitivity_pareto_plot +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab("")

new_lc_sensitivity_pareto_plot <-
  lc_sensitivity_pareto_plot +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab("                                        Fracti")

new_bc_sensitivity_pareto_plot <- 
  bc_sensitivity_pareto_plot +
  ylab("Fract on of appropriate cases recommended")

pareto_grid <- new_cvd_sensitivity_pareto_plot +
  new_bc_sensitivity_pareto_plot +
  new_lc_sensitivity_pareto_plot +
  plot_layout(ncol=1, nrow=3)

pareto_grid

ggsave(paste(figure_save_path, "appendix_figure_3.pdf", sep = ""),
       width = 12,
       height = 9)

