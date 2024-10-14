library(igraph)

calc_total_distance <- function(graph) {
  all_paths <- list()
  total_distance <- 0
  
  for (node in V(graph)$name) {
    if (node != "1") {
      paths <- all_simple_paths(graph, from="1", to=node, mode="out")
      all_paths[[node]] <- paths
      for (path in paths) {
        total_distance <- total_distance + (length(path) - 1)
      }
    }
  }
  total_distance
}

calc_node_dependency <- function(graph) {
  sum(degree(graph, mode="in"))
}

calc_min_path_length <- function(graph) {
  min(distances(graph)[1, -1])
}

calc_avg_path_length <- function(graph) {
  mean(distances(graph)[1, -1])
}

calc_num_prerequisites <- function(graph) {
  sum(degree(graph, mode="in") > 0)
}

calc_cumulative_prerequisites <- function(graph) {
  sum(eigen_centrality(graph)$vector)
}

calc_trait_isolation <- function(graph) {
  sum(degree(graph, mode="all") == 1)
}

# Betweenness Centrality
calc_betweenness_centrality <- function(graph) {
  max(betweenness(graph))
}

# Eigenvector Centrality
calc_eigenvector_centrality <- function(graph) {
  max(eigen_centrality(graph)$vector)
}

calc_graph_density <- function(graph) {
  edge_density(graph)
}

calc_all_measures<- function(data) {
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
    V(graph)$name <- as.character(1:num_nodes)
    
    # Including all measures 
    total_distance <- calc_total_distance(graph)
    node_dependency <- calc_node_dependency(graph)
    min_path_length <- calc_min_path_length(graph)
    avg_path_length <- calc_avg_path_length(graph)
    num_prerequisites <- calc_num_prerequisites(graph)
    cumulative_prerequisites <- calc_cumulative_prerequisites(graph)
    trait_isolation <- calc_trait_isolation(graph)
    betweenness_centrality <- calc_betweenness_centrality(graph)
    closeness_centrality <- calc_closeness_centrality(graph)
    eigenvector_centrality <- calc_eigenvector_centrality(graph)
    graph_density <- calc_graph_density(graph)
    
    data[data$adj_mat == adj_string, "total_distance"] <- total_distance
    data[data$adj_mat == adj_string, "node_dependency"] <- node_dependency
    data[data$adj_mat == adj_string, "min_path_length"] <- min_path_length
    data[data$adj_mat == adj_string, "avg_path_length"] <- avg_path_length
    data[data$adj_mat == adj_string, "num_prerequisites"] <- num_prerequisites
    data[data$adj_mat == adj_string, "cumulative_prerequisites"] <- cumulative_prerequisites
    data[data$adj_mat == adj_string, "trait_isolation"] <- trait_isolation
    data[data$adj_mat == adj_string, "betweenness_centrality"] <- betweenness_centrality
    data[data$adj_mat == adj_string, "eigenvector_centrality"] <- eigenvector_centrality
    data[data$adj_mat == adj_string, "graph_density"] <- graph_density
  }
  data
} 


all_data <- calc_all_measures(all_data)


model <- lm(step_transitions ~ 1, data = all_data)

# Apply forward selection
forward.selection.model <- step(model, 
                                scope = list(lower = ~1, 
                                             upper = ~total_distance + node_dependency + min_path_length + avg_path_length + num_prerequisites + cumulative_prerequisites + trait_isolation + betweenness_centrality + eigenvector_centrality + graph_density),
                                direction = "forward", 
                                trace = 0)


summary(forward.selection.model)


plotDVByDistance <- function(data, DV, DV_label) {
  average_data <- data %>%
    group_by(adj_mat, strategy, avg_path_length) %>%
    summarize(avg_DV = mean(!!sym(DV), na.rm = TRUE), .groups = 'drop')
  
  if (length(unique(data$num_nodes) > 1)) {
    num_nodes <- paste0(min(data$num_nodes), " - ", max(data$num_nodes))
  } else {
    num_nodes <- data$num_nodes[1]
  }
  
  plot <- ggplot(average_data, aes(x = avg_path_length, y = avg_DV, color = as.factor(strategy))) +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(title = paste0(DV_label, " by Average Path Length"),
         subtitle = paste0("Number of Nodes: ", num_nodes),
         x = "Average Path Length",
         y = DV_label) +
    theme_minimal() +
    scale_color_discrete(name = "Strategy") #+ 
  # scale_y_continuous(
  #   breaks = seq(0, max(average_data$avg_DV) + 1, by = 1),
  #   minor_breaks = seq(0, max(average_data$avg_DV), by = 0.5)
  # )
  
  print(plot)
  return(plot)
}

plotDVByDistance(all_data, "step_transitions", "Expected Success Rate of Learning Attempts")
