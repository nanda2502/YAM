setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## Loading libraries


library(igraph)
library(combinat)
library(ggplot2)
library(gridExtra)
library(ggraph)
library(dplyr)
library(grid)
library(cowplot)
library(lattice)


## Data preparation functions



calc_avg_path_length <- function(graph) {
  mean(distances(graph)[1, -1])
}

add_avg_path_length <- function(data) {
  unique_combinations <- data %>%
    select(num_nodes, adj_mat) %>%
    distinct()
  for (i in seq_len(nrow(unique_combinations))) {
    combination <- unique_combinations[i, ]
    adj_string <- combination[[which(colnames(unique_combinations) == "adj_mat")]]
    num_nodes <- combination[[which(colnames(unique_combinations) == "num_nodes")]]
    
    adjacency_vector <- as.numeric(unlist(strsplit(adj_string, "")))
    adjacency_matrix <- matrix(adjacency_vector, nrow = num_nodes, ncol = num_nodes, byrow = TRUE)
    graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
    V(graph )$name <- as.character(1:num_nodes)
    
    avg_path_length <- calc_avg_path_length(graph)
    
    data[data$adj_mat == adj_string, "avg_path_length"] <- avg_path_length
  }
  data
}


calculate_path_lengths_to_root <- function(graph, root) {
  total_distance <- 0
  
  for(node in V(graph)) {
    if(node != root) { # Avoid considering the root itself
      all_paths <- all_simple_paths(graph, from = node, to = root, mode = "out")
      for(path in all_paths) {
        # Sum the length of nodes in each path minus 1 to get the number of edges
        total_distance <- total_distance + (length(path) - 1)
      }
    }
  }
  
  return(total_distance)
}

add_total_distance <- function(data) {
  unique_combinations <- data %>%
    select(num_nodes, adj_mat) %>%
    distinct()
  
  for (i in seq_len(nrow(unique_combinations))) {
    combination <- unique_combinations[i, ]
    adj_string <- combination[[which(colnames(unique_combinations) == "adj_mat")]]
    num_nodes <- combination[[which(colnames(unique_combinations) == "num_nodes")]]
    
    adjacency_vector <- as.numeric(unlist(strsplit(adj_string, "")))
    adjacency_matrix <- matrix(adjacency_vector, nrow = num_nodes, ncol = num_nodes, byrow = TRUE)
    
    rev_adjacency_matrix <- t(adjacency_matrix)
    graph <- graph_from_adjacency_matrix(rev_adjacency_matrix, mode = "directed")
    V(graph)$name <- as.character(1:num_nodes)
    
    
    
    root_node <- 1
    total_distance <- calculate_path_lengths_to_root(graph, root_node)
    
    data[data$adj_mat == adj_string, "total_distance"] <- total_distance
  }
  
  return(data)
}

clean_file <- function(data) {
  sum_missing <- sum(data$adj_mat == "")
  if (sum_missing > 0) {
    print(paste0("Removing ", sum_missing, " observations with missing adjacency matrices."))
  }
  data <- data[data$adj_mat != "", ]
  
  num_nodes <- data$num_nodes[1]
  sum_nan <- sum(is.nan(data$step_payoff))
  if (sum_nan > 0) {
    print(paste0("Removing ", sum_nan, " observations that failed to converge."))
  }
  data <- data[!is.nan(data$step_payoff), ]
  sum_short <- sum(data$steps < num_nodes - 1)
  return(data)
}

read_file <- function(num_nodes) {
  data <- read.csv(gzfile(paste0("./output/expected_steps_", num_nodes, ".csv.gz")), stringsAsFactors = FALSE, colClasses = c(adj_mat = "character"))
  data$strategy <- factor(
    data$strategy,
    levels = c(
      "RandomLearning",
      "PayoffBasedLearning",
      "ProximalLearning",
      "PrestigeBasedLearning",
      "ConformityBasedLearning"
    ),
    labels = c(
      "Random",
      "Payoff",
      "Proximal",
      "Prestige",
      "Conformity"
    )
  )
  return(data)
}

average_over_replications <- function(data) {
  outcome_vars <- c("step_payoff", "step_transitions")
  grouping_vars <- c("num_nodes", "alpha", "strategy", "adj_mat", "steps")
  
  data <- data %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(across(all_of(outcome_vars), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
  
  return(data)
}

read_all <- function(numbers) {
  all_data <- lapply(numbers, function(num_nodes) {
    data <- read_file(num_nodes)
    data <- clean_file(data)
    data <- average_over_replications(data)
    data <- add_avg_path_length(data)
    print(paste0("Finished processing data for ", num_nodes, " nodes."))
    return(data)
  })
  
  all_data <- bind_rows(all_data)
  
  return(all_data)
}

average_over_lambda <- function(data) {
  outcome_vars <- c("step_payoff", "step_transitions")
  other_vars_to_retain <- c("num_nodes", "avg_path_length")
  
  lambda_values <- seq(0.1, 20, by = 0.1)
  
  subsets <- lapply(lambda_values, function(lambda) {
    data <- data %>%
      mutate(weight = dpois(steps, lambda))
    
    weighted_averages <- data %>%
      group_by(adj_mat, strategy, alpha) %>%
      summarize(
        across(all_of(outcome_vars), ~ sum(. * weight) / sum(weight)),
        across(all_of(other_vars_to_retain), first),
        lambda = lambda,
        .groups = 'drop'
      )
    
    return(weighted_averages)
  })
  
  result <- do.call(rbind, subsets)
  
  return(result)
}

add_ratios <- function(data) {
  data$rel_payoff <- rep(NA, nrow(data))
  data$rel_success <- rep(NA, nrow(data))
  for (i in 1:nrow(data)) {
    if (data$strategy[[i]] == "Random") {
      data$rel_payoff[[i]] <- 1
      data$rel_success[[i]] <- 1
      next
    }
    
    current_lambda <- data$lambda[[i]]
    current_adj_mat <- data$adj_mat[[i]]
    payoff <- data$step_payoff[[i]]
    success <- data$step_transitions[[i]]
    random_payoff <- data[data$lambda == current_lambda & data$adj_mat == current_adj_mat & data$strategy == "Random", "step_payoff"]
    random_success <- data[data$lambda == current_lambda & data$adj_mat == current_adj_mat & data$strategy == "Random", "step_transitions"]
    data$rel_payoff <- payoff / random_payoff
    data$rel_success <- success / random_success
  }
  
  return(data)
}

#this should do the same as the one above, but faster
add_ratios <- function(data) {
  data <- data %>%
    group_by(lambda, adj_mat, alpha) %>%
    mutate(random_payoff = ifelse(strategy == "Random", step_payoff, NA),
           random_success = ifelse(strategy == "Random", step_transitions, NA)) %>%
    mutate(random_payoff = max(random_payoff, na.rm = TRUE),
           random_success = max(random_success, na.rm = TRUE)) %>%
    mutate(rel_payoff = ifelse(strategy == "Random", 1, step_payoff / random_payoff),
           rel_success = ifelse(strategy == "Random", 1, step_transitions / random_success)) %>%
    ungroup() %>%
    select(-random_payoff, -random_success)  # Optional: remove the helper columns if not needed
  
  # Replace Inf and NaN from divisions by zero or NA calculations with NA
  data$rel_payoff[!is.finite(data$rel_payoff)] <- NA
  data$rel_success[!is.finite(data$rel_success)] <- NA
  
  return(data)
}


## Plotting functions



plotDVbyIV <- function(data, DV, DV_label, IV, IV_label, lambda_value, strategy_colors = NULL) {
  average_data <- data %>%
    filter(lambda == lambda_value) %>%
    mutate(
      !!sym(DV) := scales::rescale(!!sym(DV), to = c(0, 1)),
      !!sym(IV) := scales::rescale(!!sym(IV), to = c(0, 1))
      ) %>%
    group_by(adj_mat, strategy, !!sym(IV)) %>%
    summarize(avg_DV = mean(!!sym(DV), na.rm = TRUE), .groups = 'drop')
  
  if (length(unique(data$num_nodes)) > 1) {
    num_nodes <- paste0(min(data$num_nodes), " - ", max(data$num_nodes))
  } else {
    num_nodes <- data$num_nodes[1]
  }
  
  plot <- ggplot(average_data, aes_string(x = IV, y = "avg_DV", color = "strategy")) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      title = paste0(DV_label, " by ", IV_label),
      subtitle = paste0("Number of nodes: ", num_nodes,  ", Expected number of learning opportunities λ = ", lambda_value),
      x = IV_label,
      y = DV_label
    ) +
    theme_minimal() 
  
  if (!is.null(strategy_colors)) {
    plot <- plot + scale_color_manual(name = "Strategy", values = strategy_colors)
  } else {
    plot <- plot + scale_color_discrete(name = "Strategy")
  }
  
  print(plot)
  return(plot)
}

plot_graph <- function(adj_string) {
  num_nodes <- sqrt(nchar(adj_string))
  adjacency_vector <- as.numeric(unlist(strsplit(adj_string, "")))
  adjacency_matrix <- matrix(adjacency_vector, nrow = num_nodes, ncol = num_nodes, byrow = TRUE)
  graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
  plot(graph, layout = layout_as_tree(graph))
  
}

plotUnconstrained <- function(data, num_steps) {
  unconstrained_matrices <- c(
    "0100",
    "011000000",
    "0111000000000000",
    "0111100000000000000000000",
    "011111000000000000000000000000000000",
    "0111111000000000000000000000000000000000000000000",
    "0111111100000000000000000000000000000000000000000000000000000000"
  )
  
  if (length(unique(data$num_nodes) > 1)) {
    num_nodes <- paste0(min(data$num_nodes), " - ", max(data$num_nodes))
  } else {
    num_nodes <- data$num_nodes[1]
  }
  
  filtered_data <- data[data$adj_mat %in% unconstrained_matrices & data$step_payoff > 0,]
  
  payoff_plot <- ggplot(filtered_data, aes(x = lambda, y = step_payoff, color = strategy, group = strategy)) +
    geom_point(alpha = 0.1) + 
    geom_smooth(method = "lm", se = FALSE) + 
    theme_minimal() +
    labs(title = "Expected Payoff Per Step for Unconstrained Structures",
         subtitle = paste0("Number of Nodes: ", num_nodes),
         x = "Lambda",
         y = "Payoff per Step") +
    theme(legend.position = "right")
  
  print(payoff_plot)
  return(payoff_plot)
}


level_plot <- function(data, strategy, outcome, outcome_label, Y, Y_label, scale_settings,
                       lambda_resolution = 400, Y_resolution = 400) {
  
  # Filter data by strategy
  data <- data[data$strategy == strategy, ]
  
  # Define regular grid over lambda and the specified Y variable
  lambda_seq <- seq(min(data$lambda), max(data$lambda), length.out = lambda_resolution)
  Y_seq <- seq(min(data[[Y]]), max(data[[Y]]), length.out = Y_resolution)
  
  # Perform interpolation using akima's interp function
  interp_result <- akima::interp(
    x = data$lambda, 
    y = data[[Y]], 
    z = data[[outcome]],
    xo = lambda_seq, 
    yo = Y_seq, 
    linear = FALSE,
    duplicate = "mean"  # Handles duplicate (x, y) points
  )
  
  # Prepare data for plotting
  levelplot_data <- expand.grid(lambda = interp_result$x, Y_value = interp_result$y)
  levelplot_data[[outcome]] <- as.vector(interp_result$z)
  
  lambda_min <- min(interp_result$x)
  lambda_max <- max(interp_result$x)
  Y_min <- min(interp_result$y)
  Y_max <- max(interp_result$y)
  
  # Truncate edges by removing rows where lambda or Y are at their min/max
  truncated_data <- levelplot_data %>%
    filter(lambda != lambda_min & lambda != lambda_max & 
             Y_value != Y_min & Y_value != Y_max)
  
  plot <- ggplot(truncated_data, aes(x = lambda, y = Y_value)) +
    geom_tile(aes(fill = !!sym(outcome))) +
    scale_settings +
    labs(
      x = "Lambda", 
      y = Y_label,
      fill = "Outcome Value",
      title = paste("Ratio of", outcome_label, "for", strategy, "/ Random")
    ) +
    theme_minimal()
  
  # Check if a contour at z = 1 is feasible
  z_min <- min(truncated_data[[outcome]], na.rm = TRUE)
  z_max <- max(truncated_data[[outcome]], na.rm = TRUE)
  
  if (z_min <= 1 && z_max >= 1) {
    # Add contour if outcome values around 1 exist
    plot <- plot + geom_contour(aes(z = !!sym(outcome)), color = "white", linewidth = 1.5, breaks = 1, linetype = "solid")
  }
  
  return(plot)
}


get_consistent_scale <- function(all_data, outcome) {
  global_min <- min(all_data[[outcome]], na.rm = TRUE) - 2
  global_max <- max(all_data[[outcome]], na.rm = TRUE) + 2
  
  scale_fill_viridis_c(limits = c(global_min, global_max), option = "viridis", na.value = "transparent")
}

level_plot <- function(data, strategy, outcome, outcome_label, Y, Y_label, scale_settings,
                       lambda_resolution = 400, Y_resolution = 400) {
  
  library(ggplot2)
  library(mgcv)
  library(akima)
  
  # Filter data by strategy
  data <- data[data$strategy == strategy, ]
  
  # Remove any rows with missing values in relevant columns
  data <- na.omit(data[, c("lambda", Y, outcome)])
  
  # Rename Y variable to 'Y_value' to prevent conflicts
  data$Y_value <- data[[Y]]
  
  # Define regular grid over lambda and the specified Y variable for GAM smoothing
  lambda_seq <- seq(min(data$lambda), max(data$lambda), length.out = lambda_resolution)
  Y_seq <- seq(min(data$Y_value), max(data$Y_value), length.out = Y_resolution)
  
  # Create grid for predictions
  grid <- expand.grid(lambda = lambda_seq, Y_value = Y_seq)
  
  # Fit a GAM with thin plate splines for smoothing
  gam_formula <- as.formula(paste(outcome, "~ te(lambda, Y_value)"))
  gam_model <- gam(gam_formula, data = data)
  
  # Predict values on the grid (smoothed surface)
  grid[[outcome]] <- predict(gam_model, newdata = grid, type = "response")
  
  # Prepare data for plotting - Smoothed level plot
  plot <- ggplot(grid, aes(x = lambda, y = Y_value)) +
    geom_tile(aes(fill = .data[[outcome]])) +
    scale_settings +
    labs(
      x = "Lambda", 
      y = Y_label,
      fill = "Outcome Value",
      title = paste("Ratio of", outcome_label, "for", strategy, "/ Random")
    ) +
    theme_minimal()
  
  # Use akima to interpolate original data for contour lines
  if (nrow(data) > 3) {  # Ensure there are enough data points for interpolation
    interp_result <- with(data, akima::interp(
      x = lambda, y = Y_value, z = data[[outcome]],
      xo = lambda_seq, yo = Y_seq, 
      duplicate = "mean"
    ))
    
    # Check if interpolation range includes z = 1
    z_min <- min(interp_result$z, na.rm = TRUE)
    z_max <- max(interp_result$z, na.rm = TRUE)
    
    if (z_min <= 1 && z_max >= 1) {
      # Extract contour lines at z = 1 from interpolation
      contour_data <- contourLines(interp_result$x, interp_result$y, interp_result$z, levels = 1)
      if (length(contour_data) > 0) {
        contour_df <- do.call(rbind, lapply(contour_data, function(cont) {
          data.frame(lambda = cont$x, Y_value = cont$y)
        }))
        
        plot <- plot + geom_path(
          data = contour_df,
          aes(x = lambda, y = Y_value),
          color = "white", # Customize as needed
          size = 0.7
        )
      }
    }
  }
  
  return(plot)
}

s_curve <- function(x, total) {
  return(1 / (1 + exp(-15*((x/total) - 0.5))))
}

frequencies <- c(0.6, 0.3, 0.1, 0.0)
s_curve(frequencies, sum(frequencies))

## Analysis


data <- read_all(3:8)
data <- average_over_lambda(data)

#data <- add_ratios(data)

saveRDS(data, "data_processed_newconf.rds")
data <- readRDS("data_processed_newconf.rds")


plot_graph("0100000000100000000100000000100000000100000000110000000000000000")


fully_constrained <- c(
  "010001000",
  "0100001000010000",
  "0100000100000100000100000",
  "010000001000000100000010000001000000",
  "0100000001000000010000000100000001000000010000000",
  "0100000000100000000100000000100000000100000000100000000100000000"
)


payoff_harsh <- plotDVbyIV_binned(
  data[data$alpha == 0,],
  DV = "step_payoff", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Constraints on Learning",
  lambda_value = 2,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black" )
)

payoff_easy <- plotDVbyIV_binned(
  data[data$alpha == 0,],
  DV = "step_payoff", DV_label = "Performance",
  IV = "avg_path_length",  IV_label ="Constraints on Learning",
  lambda_value = 10,
  strategy_colors = c("Payoff" = "#20BF55", "Proximal" = "#FBB13C", "Prestige" = "#ED474A", "Conformity" = "#8B80F9","Random" = "black" )
)

grid.arrange(payoff_harsh, payoff_easy)

plotDVbyIV_binned <- function(data, DV, DV_label, IV, IV_label, lambda_value, strategy_colors = NULL) {
  library(dplyr)
  library(ggplot2)
  
  
  average_data <- data %>%
    filter(lambda == lambda_value) %>%
    #mutate(!!sym(IV) := scales::rescale(!!sym(IV), to = c(0, 1))) %>%    
    group_by(adj_mat, strategy) %>%
    summarize(avg_IV = mean(!!sym(IV), na.rm = TRUE),
              avg_DV = mean(!!sym(DV), na.rm = TRUE), .groups = 'drop')
  
  min_iv <- min(average_data$avg_IV, na.rm = TRUE)
  max_iv <- max(average_data$avg_IV, na.rm = TRUE)
  
  bin_centers1 <- seq(min_iv + 0.125, min_iv + 0.5, by = 0.125)
  bin_edges1 <- data.frame(
    bin_start = bin_centers1 - 0.125,
    bin_end = bin_centers1 + 0.125,
    bin_center = bin_centers1
  )
  
  # Create bins for the rest of the range (min + 0.5 to max_iv)
  bin_centers2 <- seq(min_iv + 0.625, max_iv, by = 0.5)
  bin_edges2 <- data.frame(
    bin_start = bin_centers2 - 0.25,
    bin_end = bin_centers2 + 0.25,
    bin_center = bin_centers2
  )
  
  # Combine both bin edges
  bin_edges <- rbind(bin_edges1, bin_edges2)
  
  visible_bins <- bin_edges[bin_edges$bin_center %% 0.25 == 0, ] 

  binned_data_list <- lapply(unique(average_data$strategy), function(strategy) {
    strategy_data <- average_data %>% filter(strategy == !!strategy)
    binned_data <- lapply(1:nrow(bin_edges), function(i) {
      bin_start <- bin_edges$bin_start[i]
      bin_end <- bin_edges$bin_end[i]
      bin_center <- bin_edges$bin_center[i]
      bin_data <- strategy_data %>%
        filter(avg_IV >= bin_start & avg_IV < bin_end)
      if (nrow(bin_data) > 0) {
        data.frame(
          strategy = strategy,
          bin_center = bin_center,
          binned_avg_DV = mean(bin_data$avg_DV, na.rm = TRUE)
        )
      } else {
        NULL
      }
    })
    do.call(rbind, binned_data)
  })
  
  binned_data <- do.call(rbind, binned_data_list)
  binned_data$binned_avg_DV <- scales::rescale(binned_data$binned_avg_DV, to = c(0, 1))
  
  if (length(unique(data$num_nodes)) > 1) {
    num_nodes <- paste0(min(data$num_nodes), " - ", max(data$num_nodes))
  } else {
    num_nodes <- data$num_nodes[1]
  }
  
  if (lambda_value == 2) {
    title_string <- "Fast"
  } else {
    title_string <- "Slow"
  }
  
  plot <- ggplot(binned_data, aes(x = bin_center, y = binned_avg_DV, color = strategy, shape = strategy, fill = strategy)) +
    geom_point(size = 2) +
    geom_smooth(method = "loess", span = 0.5, se = FALSE) +
    labs(
      title = title_string,
      x = IV_label,
      y = DV_label
    ) +
    theme_minimal() + 
    ylim(0, 1) +
    theme(
      axis.text.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      )
  
  if (!is.null(strategy_colors)) {
    plot <- plot + scale_fill_manual(name = "Strategy", values = strategy_colors)
    plot <- plot + scale_color_manual(name = "Strategy", values = strategy_colors)
  } else {
    plot <- plot + scale_fill_discrete(name = "Strategy")
    plot <- plot + scale_color_discrete(name = "Strategy")
  }
  
  plot <- plot + scale_shape_manual(name = "Strategy", values = 21:25) # Assign different shapes
  
  print(plot)
  return(plot)
}

plot_graph_panel <- function(df) {
  unique_adj_matrices <- df %>% 
    distinct(adj_mat, num_nodes)
  
  plot_list <- list()
  
  for (i in 1:nrow(unique_adj_matrices)) {
    adj_string <- unique_adj_matrices$adj_mat[i]
    num_nodes <- unique_adj_matrices$num_nodes[i]
    
    adjacency_vector <- as.numeric(unlist(strsplit(adj_string, "")))
    adjacency_matrix <- matrix(adjacency_vector, nrow = num_nodes, ncol = num_nodes, byrow = TRUE)
    
    g <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
    V(g)$name <- as.character(1:num_nodes)
    
    graph_plot <- ggraph(g, layout = "tree") +
      geom_edge_link(arrow = arrow(length = unit(2, 'mm'), type = "closed"), end_cap = circle(2, 'mm')) +
      geom_node_point(size = 2) +
      theme_void() +
      theme(
        plot.margin = unit(c(1,1,1,1), "pt"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.25)
      )
    
    plot_list[[i]] <- graph_plot
  }
  
  n_plots <- length(plot_list)
  n_cols <- ceiling(sqrt(n_plots))
  n_rows <- ceiling(n_plots / n_cols)
  
  n_cols <- 7
  n_rows <- 1
  
  combined_plot <- plot_grid(plotlist = plot_list, ncol = n_cols, nrow = n_rows, align = 'none')
  
  pdf_width <- n_cols * 1.5
  pdf_height <- n_rows * 1.5
  
  print(combined_plot)
  combined_plot
}


sampled_rows <- do.call(rbind, lapply(sort(unique(data$avg_path_length))[seq(1, 43, length.out = 7)], function(val) {
  data[data$avg_path_length == val, ][sample(sum(data$avg_path_length == val), 1), ]
}))


plot_graph_panel(sampled_rows)
