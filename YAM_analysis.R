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

plot_graph("0110000010000010000000000")

##### Figures #####
data_abs <- read_abs(3:8)

datap <- read_all(3:8)

data<- read_all(8)


data_abs <- readRDS("data_abs.rds")
data <- readRDS("data_merged.rds")

saveRDS(data, "data.rds")

data <- average_over_lambda(data)

data <- data %>%
  left_join(data_abs %>% select(strategy, adj_mat, slope, steps), 
            by = c("strategy", "adj_mat", "slope")) %>%
  rename(absorbing = steps.y,
         steps = steps.x) 

#data <- add_ratios(data)

data <- readRDS("data_processed_newconf.rds")

average_indegree <- function(graph) {
  mean(degree(graph, mode = "in"))
}

calculate_average_product <- function(graph) {
  distances <- distances(graph, v = 1, mode = "out")
  indegrees <- degree(graph, mode = "in")
  average_product <- mean(distances[-1] * indegrees[-1], na.rm = TRUE)
  return(average_product)
}

average_path_length_to_root <- function(graph) {
  all_paths <- list()
  
  dfs_paths <- function(node, path) {
    if (node == 1) {
      all_paths[[length(all_paths) + 1]] <<- path
      return()
    }
    for (predecessor in neighbors(graph, node, mode = "in")) {
      dfs_paths(predecessor, c(predecessor, path))
    }
  }
  
  for (node in V(graph)$name) {
    if (node != 1) {  # Start pathfinding from every node except the root
      dfs_paths(node, node)
    }
  }

  path_lengths <- sapply(all_paths, length) - 1  # Length minus 1 for number of edges
  average_length <- mean(path_lengths)
  return(average_length)
}

calculate_lock_measure <- function(graph, root = 1) {
  # Calculate the out-degree for each node
  out_degrees <- degree(graph, mode = "out")
  
  # Calculate the shortest path distances from the root node
  distances <- distances(graph, v = root, mode = "out")[1, ]
  
  # Identify unlocking nodes (those with more than 1 outgoing edge)
  unlocking_nodes <- which(out_degrees > 1)
  
  # Calculate the locking measure
  lock_measure <- sum(unlist(lapply(unlocking_nodes, function(node) {
    num_unlocked_nodes <- out_degrees[node]
    depth_from_root <- distances[node]
    num_unlocked_nodes * depth_from_root
  })))
  
  return(lock_measure)
}

unlocking_degree <- function(graph, root = 1) {
  # Calculate the out-degree for each node
  out_degrees <- degree(graph, mode = "out")
  
  return(sum(out_degrees))
}

downstream_complexity <- function(g) {
  root <- V(g)[1]
  non_root_vertices <- V(g)[-1]
  
  downstream_counts <- sapply(non_root_vertices, function(v) {
    descendants <- subcomponent(g, v, mode="out")
    length(descendants) - 1  # Exclude the node itself
  })
  
  mean(downstream_counts)
}

mean_indegree <- function(graph) {
  mean(degree(graph, mode = "in"))
}

root_outdegree <- function(graph) {
  degree(graph, 1, mode = "out")
}

avg_root_distance <- function(g) {
  vertices <- setdiff(V(g), 0)
  
  all_paths <- lapply(vertices, function(v) {
    all_simple_paths(g, from = 1, to = v)
  })

  path_lengths <- unlist(lapply(all_paths, function(paths) {
    if (length(paths) == 0) {
      return(NA)
    }
    lengths <- sapply(paths, length)
    lengths - 1
  }))
  
  mean(path_lengths, na.rm = TRUE)
}

data <- add_graph_measure(data, avg_root_distance, "total_distance")

data_1 <- add_graph_measure(data_1, root_outdegree, "root_outdegree")

data <- add_graph_measure(data, root_outdegree, "root_outdegree")

data_perf <- add_graph_measure(data_perf, mean_indegree, "mean_indegree")

data_abs <- add_graph_measure(data_abs, downstream_complexity, "downstream_complexity")

data <- add_graph_measure(data, unlocking_degree, "unlocking_degree")

data <- add_graph_measure(data, average_indegree, "avg_indegree")

data <- add_graph_measure(data, calculate_average_product, "avg_product")

data <- add_graph_measure(data, average_path_length_to_root, "avg_total_path_length")

data <- add_graph_measure(data, calculate_modularity, "modularity")

data <- add_graph_measure(data, calculate_clustering_coefficient, "clustering")

data <- add_graph_measure(data, calculate_lock_measure, "locking")

###### Figure 1C #####

payoff_fast <- plotDVbyIV_binned(
  data[data$alpha == 0,],
  DV = "step_payoff", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Constraints on Learning",
  lambda_value = 2,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black", "Perfect" = "blue" )
)

payoff_slow <- plotDVbyIV_binned(
  data[data$alpha == 0,],
  DV = "step_payoff", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Constraints on Learning",
  lambda_value = 10,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black", "Perfect" = "blue" )
)

###### Success ~ Constraints ####

success_slow <- plotDVbyIV(
  data[data$avg_path_length == 2,],
  DV = "step_transitions", DV_label = "Learning Success Rate",
  IV = "avg_path_length",  IV_label ="Constraints on Learning",
  lambda_value = 10,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black", "Perfect" = "blue" )
)

plotDVbyIV(
  get_default(data),
  DV = "step_payoff", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Average Path Length",
  lambda_value = 5,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black", "Perfect" = "blue" )
)

plotDVbyIV(
  get_default(data),
  DV = "step_transitions", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Average Path Length",
  lambda_value = 1,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black", "Perfect" = "blue" )
)

graph_ids <- data.frame(
  graph = unique(data$adj_mat),
  ID = 1:length(unique(data$adj_mat))
)

plot_graph(graph_ids$graph[graph_ids$ID == 615])

plotDVbyIV(
  get_default(data[data$strategy != "Perfect" & data$num_nodes == 8,]),
  DV = "step_payoff", DV_label = "performance",
  IV = "avg_path_length",  IV_label ="mean distance to root",
  lambda_value = 4,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black")
)

plotDVbyIV(
  get_default(data[data$strategy != "Perfect" & data$num_nodes == 8,]),
  DV = "absorbing", DV_label = "total learning time",
  IV = "avg_path_length",  IV_label ="mean distance to root",
  lambda_value = 4,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black")
)

plotDVbyIV(
  get_default(data[data$strategy != "Perfect" & data$num_nodes == 8,]),
  DV = "step_variation", DV_label = "cultural diversity",
  IV = "avg_path_length",  IV_label ="mean distance to root",
  lambda_value = 20,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black")
)

plotDVbyIV(
  get_default(data_rev[data_rev$strategy != "Perfect",]),
  DV = "step_payoff", DV_label = "performance",
  IV = "avg_path_length",  IV_label ="mean distance to root",
  lambda_value = 4,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black")
)

plotDVbyIV(
  get_default(data_alpha[data_rev$strategy != "Perfect",]),
  DV = "step_payoff", DV_label = "performance",
  IV = "avg_path_length",  IV_label ="mean distance to root",
  lambda_value = 4,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black")
)

plotDVbyIV(
  data_perf,
  DV = "step_payoff", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Average path length",
  lambda_value = NULL,
  strategy_colors = c("Perfect" = "blue" )
)

plotDVbyIV(
  default_data,
  DV = "step_transitions", DV_label = "Success rate",
  IV = "prop_learnable",  IV_label ="Average Proportion Learnable",
  lambda_value = NULL,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black", "Perfect" = "blue" )
)

plotDVbyIV(
  data_abs[data_abs$num_nodes == 8 & data_abs$strategy == "Proximal" & data_abs$slope == 2,],
  DV = "step_payoff", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Average path length",
  lambda_value = NULL
)

plotDVbyIV(
  data,
  DV = "step_payoff", DV_label = "Performance",
  IV = "root_outdegree",  IV_label ="Root Outdegree",
  lambda_value = 5
)

plotDVbyIV(
  data,
  DV = "step_payoff", DV_label = "Performance",
  IV = "scaled_outdegree",  IV_label ="Weighted Root Outdegree",
  lambda_value = 5,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black", "Perfect" = "blue" )
)

plotDVbyIV(
  data[!data$strategy %in% c("Payoff", "Conformity"),],
  DV = "absorbing", DV_label = "Total Learning Time",
  IV = "avg_path_length",  IV_label ="Average Path Length",
  lambda_value = NULL,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black", "Perfect" = "blue" )
)

plotDVbyIV(
  data,
  DV = "step_transitions", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Average path length",
  lambda_value = 1,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black", "Perfect" = "blue" )
)

plotDVbyIV_outdeg(
  data[data$strategy == "Payoff", ],
  DV = "absorbing", DV_label = "total learning time",
  IV = "avg_path_length",  IV_label ="Average path length",
  lambda_value = 1,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black", "Perfect" = "blue" )
)


plot_graph("0101000101000000000000000")

summary(lm(step_payoff ~  locking + avg_path_length + avg_indegree + strategy, data = data[data$num_nodes == 8,]))

flexplot(step_payoff ~ avg_path_length + avg_indegree, data = data[data$num_nodes == 8,], sample = 200)

cor(data[data$num_nodes == 8,c("avg_path_length", "avg_product", "avg_indegree")])


### min vs max payoff within a specified range
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

plot_slopes <- function(strategy, data) {
  plot <- ggplot(data[data$steps == 4 & data$strategy == strategy,], aes(x = avg_path_length, y = step_payoff, color = as.factor(slope), group = slope)) +
    geom_smooth(method = "loess", se = FALSE) +
    geom_smooth(data = data[data$steps == 4 & data$strategy == "Random", ],
                aes(x = avg_path_length, y = step_payoff), 
                method = "loess", se = FALSE, color = "black", linetype = "dashed") + 
    labs(x = "mean distance to root", y = "performance", color = "strength of bias") +
    ggtitle(strategy) +
    theme_minimal() + 
    ylim(min(data$step_payoff), max(data$step_payoff))
  print(plot)
  return(plot)
}

strategies <- c("Payoff", "Proximal", "Prestige", "Conformity")
plots <- vector("list", length(strategies))
for (i in 1:length(strategies)) {
  plots[[i]] <- plot_slopes(strategies[i], data[data$num_nodes == 8 & !(data$strategy == "Conformity" & data$slope == 0),])
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
cols <- c("black",  # Placeholder or a value for index 1 (won't be used)
          "red",    # Used for strat = 2
          "blue",   # Used for strat = 3
          "green",  # Used for strat = 4
          "purple"  # Used for strat = 5
)
slopes<-rbind(c(0,5,9),
              c(1,2,5),
              c(1,2,5),
              c(0,5,15))
xlims<-rbind(c(0,1),
             c(1,8),
             c(1,8),
             c(0,1))
xlabs<-c("Payoffs", "Trait difference", "Trait difference", "Frequency")
par(mfrow=c(2,2))
for (strat in 2:5){

  minX<-xlims[strat-1,1]
  maxX<-xlims[strat-1,2]

  x<-seq(from = minX, to = maxX, length.out = 7)
  plot(0, type='n', xlim=c(minX, maxX), ylim=c(0,1), xlab=xlabs[strat-1], ylab="Weight", axes=FALSE)
  axis(1)
  axis(2)

  for (k in 1:3){
    b<-slopes[strat-1,k]

    if (strat==2){ # payoff bias
      y<-x^b
    }
    if (strat==3) { # proximal learning
      y<-b^(1-x)
    }
    if (strat==4){ # prestige learning
      y<-b^(x-1)
    }
    if (strat==5){
      y<-1 / (1 + exp(-b * (x-0.5))) # conformity
    }

    y<-y/max(y)

    lines(x,y, col=cols[strat], lwd=2)

    for (i in 1:length(x)) points(x[i], y[i], pch=14+k, col=cols[strat], cex=1.5)
  }
}


s_curve <- function(x, total, offset) {
  return(1 / (1 + exp(-5*((x/total) - offset))))
}

frequencies <- c(0.5, 0.3, 0.15, 0.05)

par(mfrow = c(1, 2))
plot(s_curve(frequencies[1:4], 1, 1/length(frequencies[1:4])), type = "l", main = paste("Frequencies:", paste(frequencies, collapse = " ")), sub = "Offset = 1/length(frequencies)", xlab = "x", ylab = "y", ylim = c(0, 1))
plot(s_curve(frequencies[1:4], 1, 0.5), type = "l", sub = "Offset = 0.5", xlab = "x", ylab = "y", ylim = c(0, 1))

data_perf <- data[data$strategy == "Perfect", ]
data_perf$payoff_scaled <- scales::rescale(data_perf$step_payoff, to = c(0, 1))
plot_graph_panel(data_perf[data_perf$payoff_scaled > 0.6 & data_perf$payoff_scaled < 0.7,])
plot_graph_panel(data_perf[data_perf$payoff_scaled > 0.4 & data_perf$payoff_scaled < 0.5,])

data_1 <- read_all(99)
data_learn <- read.csv("../Cassava/results.csv", colClasses = c(adj_mat = "character"))
data_1 <- get_default(data_1)
data_learn <- get_default(data_learn)
data_abs <- merge(data_fixed, data_learn, by = c("adj_mat", "strategy", "slope"))
data_abs <- get_default(data_abs)


plotDVbyIV(
  data_abs,
  DV = "step_transitions", DV_label = "Success Rate",
  IV = "prop_learnable",  IV_label ="Proportion Learnable",
  lambda_value = NULL,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black", "Perfect" = "blue" )
)
plotDVbyIV_outdeg(
  data_abs[data_abs$strategy == "Proximal",],
  DV = "step_payoff", DV_label = "Performance",
  IV = "prop_learnable",  IV_label ="Proportion Learnable",
  lambda_value = NULL,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black", "Perfect" = "blue" )
)

graph_ids <- data.frame(
  graph = unique(data_abs$adj_mat),
  ID = 1:length(unique(data_abs$adj_mat))
)

data_abs$graph_id <- graph_ids$ID[match(data_abs$adj_mat, graph_ids$graph)]

plot_graph(data_abs$adj_mat[data_abs$graph_id == 86])


data$scaled_outdegree <- NULL
for (row in seq_len(nrow(data_1))) {
  num_nodes <- data$num_nodes[row]
  data$scaled_outdegree[row] <- data$root_outdegree[row] * node_weights[num_nodes - 2]
}
plot(
  3:8,
  tapply(
    data$scaled_outdegree[data$steps == 1] / data$root_outdegree[data$steps == 1],
    data$num_nodes[data$steps == 1],
    mean
  ),
  type = 'l',
  ylab = "Outdegree weight",
  xlab = "structure size"
)

