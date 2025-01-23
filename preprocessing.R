
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

add_graph_measure <- function(data, measure_func, measure_name) {
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
    
    measure <- measure_func(graph)
    
    data[data$adj_mat == adj_string, measure_name] <- measure
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
  data <- read.csv(paste0("./output/expected_steps_", num_nodes, ".csv"), stringsAsFactors = FALSE, colClasses = c(adj_mat = "character"))
  data$strategy <- factor(
    data$strategy,
    levels = c(
      "RandomLearning",
      "PayoffBasedLearning",
      "ProximalLearning",
      "PrestigeBasedLearning",
      "ConformityBasedLearning",
      "PerfectLearning"
    ),
    labels = c(
      "Random",
      "Payoff",
      "Proximal",
      "Prestige",
      "Conformity",
      "Perfect"
    )
  )
  return(data)
}

average_over_replications <- function(data) {
  outcome_vars <- c("step_payoff", "step_transitions", "step_variation")
  grouping_vars <- c("num_nodes", "alpha", "strategy", "adj_mat", "steps", "slope")
  
  data <- data %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(across(all_of(outcome_vars), ~ mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
    ungroup()
  
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
  outcome_vars <- c("step_payoff", "step_transitions", "step_variation")
  other_vars_to_retain <- c("num_nodes", "avg_path_length")
  
  lambda_values <- seq(0.1, 20, by = 0.1)
  
  subsets <- lapply(lambda_values, function(lambda) {
    data <- data %>%
      mutate(weight = dpois(steps, lambda))
    
    weighted_averages <- data %>%
      group_by(adj_mat, strategy, alpha, slope) %>%
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