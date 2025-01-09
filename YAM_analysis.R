##### Load Libraries #####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(igraph)
library(combinat)
library(ggplot2)
library(gridExtra)
library(ggraph)
library(dplyr)
library(grid)
library(cowplot)
library(lattice)

source("preprocess.R")
source("plotting.R")

##### Figures #####
data <- read_all(3:8)
data$expected_traits <-  data$step_transitions * data$steps
saveRDS(data, "data_raw.rds")

data <- average_over_lambda(data)
saveRDS(data, "data_lambda.rds")


data <- readRDS("data_lambda.rds")



###### Figure 1C #####

payoff_fast <- plotDVbyIV_binned(
  data[data$alpha == 0,],
  DV = "step_payoff", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Constraints on Learning",
  lambda_value = 2,
  strategy_colors = 
    c("Payoff" = "#20BF55", "Proximal" = "#FBB13C",
      "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black" )
)

payoff_slow <- plotDVbyIV_binned(
  data[data$alpha == 0,],
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "avg_path_length",
  IV_label ="Constraints on Learning",
  lambda_value = 10,
  strategy_colors = 
    c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", 
      "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black" )
)

###### Success ~ Constraints ####

success_slow <- plotDVbyIV(
  data[data$alpha == 0,],
  DV = "step_transitions", DV_label = "Learning Success Rate",
  IV = "avg_path_length",  IV_label ="Constraints on Learning",
  lambda_value = 10,
  strategy_colors = 
    c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", 
      "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black" )
)




###### Variation ~ Time ####

data_raw <- readRDS("data_raw.rds")

variation <- data_raw %>%
  filter(strategy %in% c("Payoff", "Proximal"),
         num_nodes == 8,
         step_variation < 0.005) %>%
  mutate(structure = ifelse(avg_path_length < 1.3, "Low Constraints", 
                            ifelse(avg_path_length > 3, "High Constraints", NA ))) %>%
  filter(!is.na(structure)) %>%
  mutate(structure = factor(
    structure,
    levels = c("Low Constraints", "High Constraints")
  )) %>%
  plotDVbyTime(
    DV = "step_variation", DV_label = "Cultural Variation",
    IV = "steps",  IV_label ="Time",
    strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C")
)

variation

###### Traits/Performance ~ Time ####

success <- data_raw %>%
  filter(strategy %in% c("Payoff", "Proximal"),
         num_nodes == 8) %>%
  mutate(structure = ifelse(avg_path_length < 1.3, "Low Constraints", 
                            ifelse(avg_path_length > 3, "High Constraints", NA ))) %>%
  filter(!is.na(structure)) %>%
  mutate(structure = factor(
    structure,
    levels = c("Low Constraints", "High Constraints")
  )) %>%
  plotDVbyTime(
  DV = "expected_traits", DV_label = "Number of Traits Learned",
  IV = "steps",  IV_label ="Time",
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C")
)

payoff <- data_raw %>%
  filter(strategy %in% c("Payoff", "Proximal"),
         num_nodes == 8) %>%
  mutate(structure = ifelse(avg_path_length == 1, "Low Constraints", 
                            ifelse(avg_path_length > 3, "High Constraints", NA ))) %>%
  filter(!is.na(structure)) %>%
  mutate(structure = factor(
    structure,
    levels = c("Low Constraints", "High Constraints")
  )) %>%
  plotDVbyTime(
    DV = "step_payoff", DV_label = "Performance",
    IV = "steps",  IV_label ="Time",
    strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C")
  )


data_combined <- data_raw %>%
  filter(
    strategy %in% c("Payoff", "Proximal"),
    num_nodes == 8
  ) %>%
  mutate(structure = case_when(
    avg_path_length < 1.2 ~ "Low Constraints",
    avg_path_length > 3 ~ "High Constraints",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(structure)) %>%
  mutate(structure = factor(
    structure,
    levels = c("Low Constraints", "High Constraints")
  )) %>%
  group_by(strategy, steps, structure) %>%
  summarize(
    avg_expected_traits = mean(expected_traits, na.rm = TRUE),
    avg_cum_payoff = mean(step_payoff, na.rm = TRUE),
    .groups = 'drop'
  )

min_expected_traits <- min(data_combined$avg_expected_traits, na.rm = TRUE)
max_expected_traits <- max(data_combined$avg_expected_traits, na.rm = TRUE)
range_expected_traits <- max_expected_traits - min_expected_traits

min_cum_payoff <- min(data_combined$avg_cum_payoff, na.rm = TRUE)
max_cum_payoff <- max(data_combined$avg_cum_payoff, na.rm = TRUE)
range_cum_payoff <- max_cum_payoff - min_cum_payoff

scaling_factor <- range_expected_traits / range_cum_payoff

data_combined <- data_combined %>%
  mutate(
    cum_payoff_rescaled = 
      (avg_cum_payoff - min_cum_payoff) * scaling_factor + min_expected_traits
  )

ggplot() +
  geom_line(
    data = data_combined,
    aes(
      x = steps,
      y = avg_expected_traits,
      color = strategy,
      linetype = "Number of Traits Learned"
    ),
    size = 1.2
  ) +
  geom_line(
    data = data_combined,
    aes(
      x = steps,
      y = cum_payoff_rescaled,
      color = strategy,
      linetype = "Performance"
    ),
    size = 1.2
  ) +
  scale_y_continuous(
    name = "Number of Traits Learned",
    sec.axis = sec_axis(
      trans = ~ (. - min_expected_traits) / scaling_factor + min_cum_payoff,
      name = "Performance"
    )
  ) +
  scale_color_manual(name = "Strategy",
                     values = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C")) +
  scale_linetype_manual(
    name = "Variable",
    values = c(
      "Number of Traits Learned" = "solid",
      "Performance" = "dotted"
    )
  ) +
  labs(x = "Time", title = "Performance and number of traits learned over time") +
  facet_wrap(~ structure) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



##### Extra things #####



sampled_rows <- do.call(rbind, lapply(
  sort(unique(data$avg_path_length))[seq(1, 43, length.out = 7)],
  function(val) {
    data[data$avg_path_length == val, ][sample(sum(data$avg_path_length == val), 1), ]
  }
))




plot_graph_panel(data_variation[data_variation$step_variation > 0.003 & data_variation$num_nodes == 8,])

s_curve <- function(x, total) {
  return(1 / (1 + exp(-15*((x/total) - 0.5))))
}

frequencies <- c(0.6, 0.3, 0.1, 0.0)
s_curve(frequencies, sum(frequencies))


plot_graph("0100000000100000000100000000100000000100000000110000000000000000")


fully_constrained <- c(
  "010001000",
  "0100001000010000",
  "0100000100000100000100000",
  "010000001000000100000010000001000000",
  "0100000001000000010000000100000001000000010000000",
  "0100000000100000000100000000100000000100000000100000000100000000"
)

