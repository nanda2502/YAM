plot_graph <- function(adj_mat, repertoire) {
  graph <- graph_from_adjacency_matrix(adj_mat, mode = "directed", weighted = TRUE, diag = TRUE)
  if (vcount(graph) != length(repertoire)) stop("Mismatch between nodes and repertoire length.")
  
  plot <- ggraph(graph, layout = "tree") +
    geom_node_point(aes(fill = factor(repertoire)), size = 12, color = "black", stroke = 1, shape = 21) +
    geom_edge_link(color = "black",
                   width = 3,
                   arrow = arrow(type = "open", length = unit(2, "mm")),
                   start_cap = circle(10, 'mm'),
                   end_cap = circle(10, 'mm')) +
    scale_fill_manual(values = c("white", "black")) +
    theme_void() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "transparent", color = NA),  # Make panel background transparent
      plot.background = element_rect(fill = "transparent", color = NA)   # Make plot background transparent
    )
  plot
}



plot_graph(adj_matrix, repertoires[[8]])

ggsave("graph.jpg", subplot, bg = 'transparent', width = 10, height = 10, units = "cm")

adj_matrix <- matrix(
  c(0, 1, 0, 0, 0, 1,
    0, 0, 1, 0, 1, 0,
    0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0),
  nrow = 6,
  byrow = TRUE
)


repertoires <- list(
  c(1, 0, 0, 0, 0, 0),
  c(1, 1, 0, 0, 0, 0),
  c(1, 0, 0, 0, 0, 1),
  c(1, 1, 1, 0, 0, 0),
  c(1, 1, 0, 0, 1, 0),
  c(1, 1, 0, 0, 0, 1),
  c(1, 1, 1, 1, 0, 0),
  c(1, 1, 1, 0, 1, 0),
  c(1, 1, 1, 0, 0, 1),
  c(1, 1, 0, 0, 1, 1),
  c(1, 1, 1, 1, 1, 0),
  c(1, 1, 1, 1, 0, 1),
  c(1, 1, 1, 0, 1, 1),
  c(1, 1, 1, 1, 1, 1)
)

transition_mat <- matrix(
  c(
    0.39, 0.48, 0.13, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.01, 0.00, 0.02, 0.71, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.45, 0.00, 0.00, 0.55, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.01, 0.73, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.02, 0.00, 0.00, 0.08, 0.00, 0.89, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.01, 0.00, 0.00, 0.03, 0.96, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.74, 0.26, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.02, 0.00, 0.98, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.01, 0.99, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.21, 0.00, 0.00, 0.79, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00
  ),
  nrow = 14,
  byrow = TRUE
)

plot_transition_mat <- function(adj_matrix) {
  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = TRUE, diag = TRUE)
  
  edge_weights <- E(graph)$weight
  
  ggraph(graph, layout = "tree") +
    geom_edge_link(
      aes(width = edge_weights),
      color = "forestgreen",
      alpha = 1,
      arrow = arrow(type = "closed", length = unit(2, "mm")),
      start_cap = circle(5, 'mm'),
      end_cap = circle(5, 'mm')
    ) +
    geom_edge_loop(
      aes(width = edge_weights),
      color = "red",
      arrow = arrow(type = "closed", length = unit(2, "mm")),
      alpha = 0.5,
      start_cap = circle(5, 'mm'),
      end_cap = circle(5, 'mm')
    ) +
    geom_node_point(size = 1, fill = "white", color = "black", stroke = 1, shape = 21) +
    scale_edge_width(range = c(0.1, 2)) +
    theme_void() +
    theme(legend.position = "none")
}

plot_transition_mat(transition_mat)




