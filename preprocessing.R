
calc_avg_path_length <- function(graph) {
  mean(distances(graph)[1, -1])
}

add_avg_path_length <- function(data) {
  unique_combinations <- data %>%
    dplyr::select(adj_mat) %>%
    distinct()
  for (i in seq_len(nrow(unique_combinations))) {
    combination <- unique_combinations[i, ]
    adj_string <- combination[[which(colnames(unique_combinations) == "adj_mat")]]
   
    graph <- string_to_igraph(adj_string)
    V(graph )$name <- as.character(1:vcount(graph))
    
    avg_path_length <- calc_avg_path_length(graph)
    
    data[data$adj_mat == adj_string, "avg_path_length"] <- avg_path_length
  }
  data
}

add_weighted_r <- function(data) {
  unique_combinations <- data %>%
    dplyr::select(adj_mat, edge_weight) %>%
    distinct()
  
  data$weighted_r <- NULL
  
  results <- data.frame(
    adj_mat = character(),
    edge_weight = numeric(),
    weighted_r = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(unique_combinations))) {
    combination <- unique_combinations[i, ]
    adj_string <- combination[[which(colnames(unique_combinations) == "adj_mat")]]
    g <- string_to_igraph(adj_string)
    V(g)$name <- as.character(1:vcount(g))
    
    edge_weight <- combination[[which(colnames(unique_combinations) == "edge_weight")]]
    
    E(g)$weight <- edge_weight
    
    nodes <- setdiff(V(g), V(g)[1])
    
    sum_weights <- sapply(nodes, function(v) {
      in_edges <- incident(g, v, mode = "in")
      
      edge_sources <- ends(g, in_edges)[, 1]  # Get the source nodes of these edges
      
      non_root_edges <- in_edges[edge_sources != 1]
      
      if (length(non_root_edges) == 0) {
        return(0)
      }
      
      edge_weights <- E(g)[non_root_edges]$weight
      
      return(sum(edge_weights))
    })
    

    results <- rbind(results, data.frame(
      adj_mat = adj_string,
      edge_weight = edge_weight,
      weighted_r = mean(sum_weights)
    ))
  }
  data <- data %>%
    left_join(results, by = c("adj_mat", "edge_weight"))
  data
}

add_graph_measure <- function(data, measure_func, measure_name) {
  unique_combinations <- data %>%
    dplyr::select(adj_mat) %>%
    distinct()
  for (i in seq_len(nrow(unique_combinations))) {
    combination <- unique_combinations[i, ]
    adj_string <- combination[[which(colnames(unique_combinations) == "adj_mat")]]
    num_nodes <- nchar(adj_string)

    graph <- string_to_igraph(adj_string)
    V(graph )$name <- as.character(1:vcount(graph))
    
    measure <- measure_func(graph)
    
    data[data$adj_mat == adj_string, measure_name] <- measure
  }
  data
}


mean_prereq <- function(g) {
  r <- if ("name" %in% vertex_attr_names(g)) V(g)[V(g)$name == "1"] else V(g)[1]
  s <- if ("name" %in% vertex_attr_names(g)) V(g)[V(g)$name != "1"] else V(g)[-1]
  mean(sapply(s, function(v) length(setdiff(subcomponent(g, v, mode = "in"), c(v, r)))))
}

mean_prereq_tc <- function(g) {
  #mean indegree of all nodes except the root
  g <- induced_subgraph(g, V(g)[-1])
  r_vals <- igraph::strength(g, mode = "in")
  mean(r_vals)
}

calculate_path_lengths_to_root <- function(graph, root = 0) {
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
    dplyr::select(num_nodes, adj_mat) %>%
    distinct()
  
  for (i in seq_len(nrow(unique_combinations))) {
    combination <- unique_combinations[i, ]
    adj_string <- combination[[which(colnames(unique_combinations) == "adj_mat")]]

    graph <- string_to_igraph(adj_string)
    V(graph )$name <- as.character(1:vcount(graph))
    
    
    
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
  return(data)
}

average_over_replications <- function(data) {
  outcome_vars <- c("step_payoff", "step_transitions", "step_variation")
  grouping_vars <- c("num_nodes", "alpha", "strategy", "adj_mat", "steps", "slope", "distribution", "lambda")
  
  data <- data %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(across(all_of(outcome_vars), ~ mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
    ungroup()
  
  return(data)
}

add_expected_traits <- function(data) {
  data %>%
    group_by(adj_mat, strategy, slope, distribution) %>%
    arrange(steps, .by_group = TRUE) %>%
    mutate(expected_traits = cumsum(step_transitions)) %>%
    ungroup()
}

read_abs <- function(numbers) {
  all_data <- lapply(numbers, function(num_nodes) {
    data <- read_file(num_nodes)
    data <- clean_file(data)
    data <- add_avg_path_length(data)
    print(paste0("Finished processing data for ", num_nodes, " nodes."))
    return(data)
  })
  
  all_data <- bind_rows(all_data)
  
  return(all_data)
}

read_sim <- function(numbers) {
  all_data <- lapply(numbers, function(num_nodes) {
    data <- read_file(num_nodes)
    data <- clean_file(data)
    data <- add_avg_path_length(data)
    data <- add_graph_measure(data, mean_prereq, "mean_prereq")
    print(paste0("Finished processing data for ", num_nodes, " nodes."))
    return(data)
  })
  
  all_data <- bind_rows(all_data)
  
  return(all_data)
} 

read_all <- function(numbers) {
  all_data <- lapply(numbers, function(num_nodes) {
    data <- read_file(num_nodes)
    data <- clean_file(data)
    data <- add_avg_path_length(data)
    data <- add_expected_traits(data)
    data_unweighted <- data[data$edge_weight == 1,]
    data_weighted <- data[data$edge_weight < 1,]
    rm(data)
    data_unweighted <- add_graph_measure(data_unweighted, mean_prereq, "mean_prereq")
    data_weighted <- add_graph_measure(data_weighted, mean_prereq_tc, "mean_prereq")
    data <- rbind(data_unweighted, data_weighted)
    print(paste0("Finished processing data for ", num_nodes, " nodes."))
    return(data)
  })
  
  all_data <- bind_rows(all_data)
  
  return(all_data)
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

get_default <- function(data) {
  default_slopes <- list(
    "Payoff" = 2.0,
    "Proximal" = 2.0,
    "Prestige" = 2.0,
    "Conformity" = 2.0,
    "Random" = 0.0,
    "Perfect" = 0.0
  )
  
  data %>% 
    filter(
      slope == sapply(as.character(strategy), function(x) default_slopes[[x]]),
      distribution == "Learnability",
      payoffdist == 0,
      alpha == 0,
      edge_weight == 1,
      lambda == 0
           )
  
}


get_default_slopes <- function(data) {
  default_slopes <- list(
    "Payoff" = 2.0,
    "Proximal" = 2.0,
    "Prestige" = 2.0,
    "Conformity" = 2.0,
    "Random" = 0.0,
    "Perfect" = 0.0
  )
  
  data %>% 
    filter(
      slope == sapply(as.character(strategy), function(x) default_slopes[[x]])
    )
  
}

string_to_igraph <- function(adj_string) {
  # Check if the string contains commas (comma-separated format)
  if (grepl(",", adj_string)) {
    # Parse comma-separated format
    values <- as.numeric(unlist(strsplit(adj_string, ",")))
    n <- sqrt(length(values))
    if (n != floor(n)) {
      stop("Number of values does not form a square matrix")
    }
    n <- as.integer(n)
    
    adj_matrix <- matrix(values, nrow = n, byrow = TRUE)
    weighted <- TRUE
  } else {
    # Get the length of the string
    str_length <- nchar(adj_string)
    
    # Calculate the dimension of the square matrix
    n <- sqrt(str_length)
    if (n != floor(n)) {
      stop("Input string length is not a perfect square")
    }
    n <- as.integer(n)
    
    # Initialize adjacency matrix
    adj_matrix <- matrix(0, nrow = n, ncol = n)
    
    # Parse the string and fill the adjacency matrix
    str_chars <- strsplit(adj_string, "")[[1]]
    
    # Check if it's a binary string (contains only 0s and 1s)
    is_binary <- all(str_chars %in% c("0", "1"))
    
    for (i in 1:n) {
      for (j in 1:n) {
        index <- (i-1)*n + j
        if (is_binary) {
          # Binary interpretation (0 or 1)
          adj_matrix[i, j] <- as.numeric(str_chars[index])
        } else {
          # Convert single digit to a weight (0-9 → 0.0-0.9)
          adj_matrix[i, j] <- as.numeric(str_chars[index]) / 10
        }
      }
    }
    
    weighted <- !is_binary
  }
  
  # Create igraph object from adjacency matrix
  g <- graph_from_adjacency_matrix(
    adjmatrix = adj_matrix,
    mode = "directed",
    weighted = weighted
  )
  
  # Add a graph attribute indicating whether it's weighted
  g$is_weighted <- weighted
  
  return(g)
}
add_downstream_edges <- function(adj_matrix) {
  # Make sure the input is a matrix
  if (!is.matrix(adj_matrix)) {
    stop("Input must be a matrix")
  }
  
  # Check if the matrix is square
  n <- nrow(adj_matrix)
  if (n != ncol(adj_matrix)) {
    stop("Adjacency matrix must be square")
  }
  
  # Convert adjacency matrix to igraph object
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = NULL)
  
  # Create a reachability matrix using igraph's distances function
  # For each pair of vertices, calculate if there's a path between them
  dist_matrix <- distances(g, mode = "out")
  
  # Convert the distance matrix to a binary reachability matrix
  # If distance is finite (not Inf), there's a path from i to j
  reach_matrix <- matrix(0, nrow = n, ncol = n)
  reach_matrix[is.finite(dist_matrix)] <- 1
  
  # Ensure the diagonal is 0 (no self-loops in a DAG)
  diag(reach_matrix) <- 0
  
  return(reach_matrix)
}

process_adjacency_csv <- function(input_file, output_file = "adj_out.csv", chunk_size = 100) {
  # Read the CSV file with character class
  adj_strings <- read.csv(input_file, header = FALSE, colClasses = "character", 
                          stringsAsFactors = FALSE, quote = "")
  
  # Convert to a vector for easier processing
  adj_strings <- adj_strings$V1
  
  # Set up parallel processing
  plan(multisession, workers = availableCores() - 1)
  
  # Process in chunks for better memory management
  num_strings <- length(adj_strings)
  chunks <- split(adj_strings, ceiling(seq_along(adj_strings) / chunk_size))
  
  results <- data.frame(processed = character(), stringsAsFactors = FALSE)
  
  for (chunk in chunks) {
    # Process the chunk in parallel using furrr
    chunk_results <- future_map(chunk, function(adj_string) {
      tryCatch({
        # Parse the adjacency string into a matrix
        n <- sqrt(nchar(adj_string))
        
        # Ensure n is an integer
        if (n != floor(n)) {
          stop("String length is not a perfect square: ", nchar(adj_string))
        }
        
        n <- as.integer(n)
        adjacency_vector <- as.numeric(unlist(strsplit(adj_string, "")))
        adjacency_matrix <- matrix(adjacency_vector, nrow = n, ncol = n, byrow = TRUE)
        
        # Apply the add_downstream_edges function
        result_matrix <- add_downstream_edges(adjacency_matrix)
        
        # Convert the result matrix back to a string - transposing first to maintain correct order
        result_string <- paste(as.vector(t(result_matrix)), collapse = "")
        
        return(list(processed = result_string))
      }, error = function(e) {
        warning("Error processing adjacency string: ", e$message)
        return(list(processed = NA))
      })
    }, .options = furrr_options(seed = TRUE))
    
    # Convert list of results to data frame
    chunk_df <- do.call(rbind, lapply(chunk_results, function(res) {
      data.frame(processed = res$processed, stringsAsFactors = FALSE)
    }))
    
    # Append to results
    results <- rbind(results, chunk_df)
  }
  
  # Write the results to a CSV file without header or quotes
  write.table(results, file = output_file, sep = ",", row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
}
