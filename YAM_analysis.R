##### Load Libraries #####

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

source("preprocessing.R")
source("plotting.R")


##### Figures #####
data <- read_all(3:8)

data <- readRDS("data.rds")

data_abs <- readRDS("data_absorbing.rds")

default_slopes <- list(
 "Payoff" = 5.0,
 "Proximal" = 2.0,
 "Prestige" = 2.0,
 "Conformity" = 5.0,
 "Random" = 0.0
)

data <- data %>% 
  filter(slope == sapply(as.character(strategy), function(x) default_slopes[[x]]))

saveRDS(data, "data_merged.rds")

data <- average_over_lambda(data)

#data <- add_ratios(data)
data <- readRDS("data_processed_newconf.rds")

average_indegree <- function(graph) {
  mean(degree(graph, mode = "in"))
}

data <- add_graph_measure(data, average_indegree, "avg_indegree")

data <- data %>%
  left_join(data_abs %>% select(strategy, adj_mat, steps), by = c("strategy", "adj_mat")) %>%
  rename(absorption = steps.y, steps = steps.x)

###### Figure 1C #####

payoff_fast <- plotDVbyIV_binned(
  data[data$alpha == 0,],
  DV = "step_payoff", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Constraints on Learning",
  lambda_value = 2,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black" )
)

payoff_slow <- plotDVbyIV_binned(
  data[data$alpha == 0,],
  DV = "step_payoff", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Constraints on Learning",
  lambda_value = 10,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black" )
)

###### Success ~ Constraints ####

success_slow <- plotDVbyIV(
  data[data$avg_path_length == 2,],
  DV = "step_transitions", DV_label = "Learning Success Rate",
  IV = "avg_path_length",  IV_label ="Constraints on Learning",
  lambda_value = 10,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black" )
)

plotDVbyIV(
  data,
  DV = "step_transitions", DV_label = "Success Rate",
  IV = "avg_path_length",  IV_label ="Average Path Length",
  lambda_value = 5,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black" )
)

plotDVbyIVnofilter(
  data_abs,#[data_abs$strategy != "Payoff",],
  DV = "step_transitions", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Avg Path Length",
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black" )
)

plot_graph_panel(data[data$avg_path_length == 2 & data$num_nodes == 8 & data$strategy == "Proximal" & data$step_payoff > 4, ])

data2 <- data$step_payoff[data$avg_path_length == 2 & data$steps == 5 & data$num_nodes == 8]

max_perf <- max(data$step_payoff[data$avg_path_length == 2 & data$strategy == "Proximal" & data$steps == 5 & data$num_nodes == 8])
min_perf <- min(data$step_payoff[data$avg_path_length == 2 & data$strategy == "Proximal" & data$steps == 5 & data$num_nodes == 8])

par(mfrow = c(1, 2))
data$adj_mat[data$avg_path_length == 2 & data$strategy == "Proximal" & data$steps == 5 & data$step_payoff == max_perf] %>%
  plot_graph()

data$adj_mat[data$avg_path_length == 2 & data$strategy == "Proximal" & data$steps == 5 & data$step_payoff == min_perf] %>%
  plot_graph()

plot_graphs <- function(data, strategy) {
  max_perf <- max(data$step_payoff[data$avg_path_length == 2 & data$strategy == strategy & data$steps == 5 & data$num_nodes == 8])
  min_perf <- min(data$step_payoff[data$avg_path_length == 2 & data$strategy == strategy & data$steps == 5 & data$num_nodes == 8])
  
  par(mfrow = c(1, 2))
  data$adj_mat[data$avg_path_length == 2 & data$strategy == strategy & data$steps == 5 & data$step_payoff == max_perf] %>%
    plot_graph()
  mtext("Best", side = 3, line = 0.5, cex = 1.2, at = -0.5)
  data$adj_mat[data$avg_path_length == 2 & data$strategy == strategy & data$steps == 5 & data$step_payoff == min_perf] %>%
    plot_graph()
  mtext("Worst", side = 3, line = 0.5, cex = 1.2, at = -0.5)
}

plot_graphs(data, "Random")

###### Variation ~ Time ####

data_variation <- readRDS("data_variation.rds")

variation <- data_variation %>%
  filter(strategy %in% c("Payoff", "Proximal"),
         num_nodes == 8,
         step_variation < 0.005) %>%
  mutate(structure = ifelse(avg_path_length < 1.3, "Low Constraints", 
                            ifelse(avg_path_length > 3, "High Constraints", NA ))) %>%
  filter(!is.na(structure)) %>%
  mutate(structure = factor(structure, levels = c("Low Constraints", "High Constraints"))) %>%
  plotDVbyTime(
    DV = "step_variation", DV_label = "Cultural Variation",
    IV = "steps",  IV_label ="Time",
    strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C")
  )

variation

###### Traits/Performance ~ Time ####

data_traits <- readRDS("data_traits.rds")

success <- data_traits %>%
  filter(strategy %in% c("Payoff", "Proximal"),
         num_nodes == 8) %>%
  mutate(structure = ifelse(avg_path_length < 1.3, "Low Constraints", 
                            ifelse(avg_path_length > 3, "High Constraints", NA ))) %>%
  filter(!is.na(structure)) %>%
  mutate(structure = factor(structure, levels = c("Low Constraints", "High Constraints"))) %>%
  plotDVbyTime(
    DV = "expected_traits", DV_label = "Number of Traits Learned",
    IV = "steps",  IV_label ="Time",
    strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C")
  )

payoff <- data_traits %>%
  filter(strategy %in% c("Payoff", "Proximal"),
         num_nodes == 8) %>%
  mutate(structure = ifelse(avg_path_length == 1, "Low Constraints", 
                            ifelse(avg_path_length > 3, "High Constraints", NA ))) %>%
  filter(!is.na(structure)) %>%
  mutate(structure = factor(structure, levels = c("Low Constraints", "High Constraints"))) %>%
  plotDVbyTime(
    DV = "step_payoff", DV_label = "Performance",
    DV = "step_payoff", DV_label = "Performance",
    IV = "steps",  IV_label ="Time",
    strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C")
  )


data_combined <- data_traits %>%
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
  mutate(structure = factor(structure, levels = c("Low Constraints", "High Constraints"))) %>%
  group_by(strategy, steps, structure) %>%
  summarize(
    avg_expected_traits = mean(expected_traits, na.rm = TRUE),
    avg_cum_payoff = mean(step_payoff, na.rm = TRUE),
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
    cum_payoff_rescaled = (avg_cum_payoff - min_cum_payoff) * scaling_factor + min_expected_traits
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
  scale_color_manual(name = "Strategy",
                     values = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C")) +
  scale_linetype_manual(
    name = "Variable",
    values = c(
      "Number of Traits Learned" = "solid",
      "Performance" = "dotted"
    ),
    values = c(
      "Number of Traits Learned" = "solid",
      "Performance" = "dotted"
    )
  ) +
  labs(x = "Time", title = "Performance and number of traits learned over time") +
  labs(x = "Time", title = "Performance and number of traits learned over time") +
  facet_wrap(~ structure) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



##### Extra things #####

## varying the slopes 

plot_slopes <- function(strategy) {
  plot <- ggplot(data[data$lambda == 5 & data$strategy == strategy,], aes(x = avg_path_length, y = step_payoff, color = as.factor(slope), group = slope)) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(x = "Average Path Length", y = "Performance", color = "Slope") +
    ggtitle(strategy) +
    theme_minimal()
  print(plot)
  return(plot)
}

strategies <- c("Payoff", "Proximal", "Prestige", "Conformity")
plots <- vector("list", length(strategies))
for (i in 1:length(strategies)) {
  plots[[i]] <- plot_slopes(strategies[i])
}

grid.arrange(grobs = plots, ncol = 2)


sampled_rows <- do.call(rbind, lapply(sort(unique(data$avg_path_length))[seq(1, 43, length.out = 7)], function(val) {
  data[data$avg_path_length == val, ][sample(sum(data$avg_path_length == val), 1), ]
}))




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