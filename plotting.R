plotDVbyIV <- function(data, DV, DV_label, IV, IV_label, lambda_value, strategy_colors = NULL) {
  if (!is.null(lambda_value)) {
    data <- data %>% filter(steps == lambda_value)
  }
  
  graph_ids <- data.frame(
    graph = unique(data$adj_mat),
    ID = 1:length(unique(data$adj_mat))
  )
  
  average_data <- data %>%
    # mutate(
    #   across(all_of(DV), ~scales::rescale(.x, to = c(0, 1)))
    # ) %>%
    group_by(adj_mat, strategy, !!sym(IV)) %>%
    summarize(avg_DV = mean(!!sym(DV), na.rm = TRUE), .groups = 'drop')
  

  plot <- ggplot(average_data, aes_string(x = IV, y = "avg_DV", color = "strategy", fill = "strategy")) +
    #geom_point(alpha = 0.2) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      x = IV_label,
      y = DV_label
    ) +
    theme_minimal() 
  #+ geom_text(aes(label = ID), hjust = -0.2)
  
  if (!is.null(strategy_colors)) {
    plot <- plot + scale_color_manual(name = "Strategy", values = strategy_colors) + scale_fill_manual(name = "Strategy", values = strategy_colors)
  } else {
    plot <- plot + scale_color_discrete(name = "Strategy") + scale_fill_discrete(name = "Strategy")
  }
  
  print(plot)
  return(plot)
}

plotDVbyIVTwoDatasets <- function(data1, data2, DV, DV_label, IV, IV_label, 
                                  lambda_value = NULL, 
                                  strategy_colors = NULL,
                                  dataset_labels = c("Markov Chain", "Simulations"),
                                  dataset1_linetype = "solid",
                                  dataset2_linetype = "dashed") {
  
  # Validate linetypes
  valid_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  if (!dataset1_linetype %in% valid_linetypes) {
    warning("Invalid dataset1_linetype. Using 'solid' instead.")
    dataset1_linetype <- "solid"
  }
  if (!dataset2_linetype %in% valid_linetypes) {
    warning("Invalid dataset2_linetype. Using 'dashed' instead.")
    dataset2_linetype <- "dashed"
  }
  
  # Extract only the columns we need for each dataset
  required_cols <- c(IV, DV, "strategy", "adj_mat")
  if (!is.null(lambda_value)) required_cols <- c(required_cols, "steps")
  
  # Ensure all required columns exist in each dataset
  missing_cols1 <- setdiff(required_cols, names(data1))
  missing_cols2 <- setdiff(required_cols, names(data2))
  
  if (length(missing_cols1) > 0) {
    stop("Missing required columns in data1: ", paste(missing_cols1, collapse = ", "))
  }
  if (length(missing_cols2) > 0) {
    stop("Missing required columns in data2: ", paste(missing_cols2, collapse = ", "))
  }
  
  # Select only required columns and add dataset identifier
  data1_subset <- data1[, required_cols, drop = FALSE]
  data1_subset$dataset <- dataset_labels[1]
  
  data2_subset <- data2[, required_cols, drop = FALSE]
  data2_subset$dataset <- dataset_labels[2]
  
  # Apply lambda filtering if specified
  if (!is.null(lambda_value)) {
    data1_subset <- data1_subset %>% filter(steps == lambda_value)
    data2_subset <- data2_subset %>% filter(steps == lambda_value)
  }
  
  # Combine the datasets
  combined_data <- rbind(data1_subset, data2_subset)
  
  # Process the combined data
  average_data <- combined_data %>%
    group_by(adj_mat, strategy, dataset, !!sym(IV)) %>%
    summarize(avg_DV = mean(!!sym(DV), na.rm = TRUE), .groups = 'drop')
  
  # Ensure no NA values
  average_data <- average_data[!is.na(average_data[[IV]]) & !is.na(average_data$avg_DV), ]
  
  # Convert dataset to factor to ensure consistent ordering
  average_data$dataset <- factor(average_data$dataset, levels = dataset_labels)
  
  # Create the plot
  plot <- ggplot(average_data, aes_string(x = IV, y = "avg_DV", color = "strategy", fill = "strategy")) +
    geom_smooth(aes(linetype = dataset), method = "loess", se = FALSE) +
    labs(
      x = IV_label,
      y = DV_label,
      linetype = "Dataset"
    ) +
    theme_minimal() +
    scale_linetype_manual(
      values = c(dataset1_linetype, dataset2_linetype),
      labels = dataset_labels
    )
  
  # Apply strategy colors if provided
  if (!is.null(strategy_colors)) {
    plot <- plot + 
      scale_color_manual(name = "Strategy", values = strategy_colors) + 
      scale_fill_manual(name = "Strategy", values = strategy_colors)
  } else {
    plot <- plot + 
      scale_color_discrete(name = "Strategy") + 
      scale_fill_discrete(name = "Strategy")
  }
  
  # Add a more informative legend
  plot <- plot + 
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )
  
  print(plot)
  return(plot)
}

plotDVbyIVBinned <- function(data, DV, DV_label, IV, IV_label,
                             lambda_ratio = NULL,
                             num_bins = 7,
                             DV_trans = identity,
                             log_scale_y = FALSE,
                             DV_scale = NULL,
                             # Optional bins and x-axis positions that will be generated dynamically
                             bins = NULL,
                             xposs = NULL,
                             title = NULL,
                             show_plot = FALSE,
                             # Parameters for confidence intervals
                             show_ci = FALSE,
                             conf_level = 0.95,
                             show_ylab = TRUE,
                             legend_position = "none") {
  
  # Hardcoded strategy order and associated slope values.
  # (Random = 0, all others = 2)
  strategyOrder <- c("Random", "Payoff", "Proximal", "Prestige", "Conformity")
  selectedSlopes <- c(0, 2, 2, 2, 2)
  # Optionally filter by lambda_ratio. This determines which time step to choose
  # so that the time step is proportional to the size of the tree.
  if (!is.null(lambda_ratio)) {
    data <- data %>% filter(steps == floor(lambda_ratio * num_nodes))
  }
  
  if (!is.null(DV_scale)) {
    data[[DV]] <- data[[DV]] / (DV_scale * data$num_nodes)
  }
  
  if (is.null(bins)) {
    min_iv <- min(data[[IV]], na.rm = TRUE)
    max_iv <- max(data[[IV]], na.rm = TRUE)
    
    # Small epsilon for narrow min/max bins
    epsilon <- 0.001
    
    # Create equidistant x-positions
    x_positions <- seq(min_iv, max_iv, length.out = num_bins)
    
    # Create bin boundaries with narrow min/max bins
    bins <- c(min_iv - epsilon, min_iv + epsilon)
    
    # Add middle bin boundaries
    if (num_bins > 2) {
      for (i in 2:(num_bins-1)) {
        middle_point <- (x_positions[i] + x_positions[i+1]) / 2
        bins <- c(bins, middle_point)
      }
    }
    
    # Add narrow max bin
    bins <- c(bins, max_iv - epsilon, max_iv + epsilon)
    
    # Set x-positions
    xposs <- x_positions
  }
  # Initialize a data frame to store aggregated results.
  agg_data <- data.frame()
  
  # Loop over each strategy.
  for (i in seq_along(strategyOrder)) {
    strat <- strategyOrder[i]
    slope_val <- selectedSlopes[i]
    
    # Subset the data for the current strategy and slope.
    strat_data <- data %>% 
      filter(strategy == strat, slope == slope_val)
    
    # Process each bin.
    for (j in seq_len(length(bins) - 1)) {
      bin_lower <- bins[j]
      bin_upper <- bins[j + 1]
      
      # Subset rows where the IV falls inside the current bin interval.
      bin_data <- strat_data %>% 
        filter(.data[[IV]] > bin_lower, .data[[IV]] <= bin_upper)
      
      if (nrow(bin_data) == 0) next
      
      # For each unique network in the bin, choose the DV value 
      # from the middle observation.
      networks <- unique(bin_data$adj_mat)
      mid_vals <- c()
      for (net in networks) {
        net_data <- bin_data %>% filter(adj_mat == net)
        if (nrow(net_data) == 0) next
        mid_index <- ceiling(nrow(net_data) / 2)
        mid_vals <- c(mid_vals, DV_trans(net_data[[DV]][mid_index]))
      }
      
      if (length(mid_vals) > 0) {
        avg_val <- mean(mid_vals, na.rm = TRUE)
        
        # Calculate confidence intervals
        if (length(mid_vals) >= 2) {
          # Use t-distribution for confidence interval
          t_val <- qt((1 + conf_level) / 2, df = length(mid_vals) - 1)
          sd_val <- sd(mid_vals, na.rm = TRUE)
          se_val <- sd_val / sqrt(length(mid_vals))
          ci_lower <- avg_val - t_val * se_val
          ci_upper <- avg_val + t_val * se_val
        } else {
          # If only one observation, can't calculate CI
          ci_lower <- NA
          ci_upper <- NA
        }
        
        # Use the provided x position for this bin.
        agg_data <- rbind(agg_data, data.frame(
          strategy = strat,
          bin = j,
          x = xposs[j],
          DV_val = avg_val,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          n = length(mid_vals)
        ))
      }
    }
  }
  
  # For legend ordering, we want the legend to show:
  # "Payoff", "Proximal", "Prestige", "Conformity", "Random"
  agg_data$strategy <- factor(agg_data$strategy,
                              levels = c("Payoff", "Proximal", "Prestige", "Conformity", "Random"))
  
  # Define manual mappings:
  col_map <- c(
    "Random" = "grey30",
    "Payoff" = "#006328",
    "Proximal" = "#ff8954",
    "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
    "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
  )
  shape_map <- c(
    "Payoff" = 16,
    "Proximal" = 17,
    "Prestige" = 18,
    "Conformity" = 15,
    "Random" = NA  # No shape for Random
  )
  lty_map <- c(
    "Random" = 2,
    "Payoff" = 1,
    "Proximal" = 1,
    "Prestige" = 1,
    "Conformity" = 1
  )
  
  shape_size_map <- c(
    "Random" = 3,
    "Payoff" = 3,  
    "Proximal" = 3, 
    "Prestige" = 3,
    "Conformity" = 3
  )
  
  # Create the ggplot.
  p <- ggplot(agg_data, aes(x = x, y = DV_val,
                            group = strategy,
                            color = strategy,
                            shape = strategy,
                            linetype = strategy)) +
    # Use the size mapping for line thickness
    geom_line(size = 1)
  
  # Create a subset of data for non-Random strategies (for points)
  non_random_data <- subset(agg_data, strategy != "Random")
  
  # Add white backing points first (larger size for proper coverage)
  p <- p + geom_point(
    data = non_random_data,
    aes(size = strategy),
    color = "white",
    fill = "white",
    stroke = 2  # White border/stroke for better coverage
  )
  
  # Then add the actual colored points on top
  p <- p + geom_point(
    data = non_random_data,
    aes(size = strategy, color = strategy, shape = strategy)
  )
  
  # Apply themes and scales
  p <- p + theme_classic() +
    scale_color_manual(
      values = col_map,
      breaks = strategyOrder
    ) +
    scale_shape_manual(
      values = shape_map,
      breaks = setdiff(strategyOrder, "Random")  # Exclude Random from shape legend
    ) +
    scale_linetype_manual(
      values = lty_map,
      breaks = strategyOrder
    ) +
    scale_size_manual(
      values = shape_size_map,
      breaks = strategyOrder,
      guide = "none"  # Hide the size legend
    ) +
    labs(x = IV_label, 
         y = DV_label
    ) +
    scale_x_continuous(
      limits = c(0, max(agg_data$x) * 1.05),  # Add some padding on the right
      expand = expansion(mult = c(0.05, 0.05))  # Add consistent padding on both sides
    ) +
    theme(
      # Increase the thickness of the x and y axis lines.
      axis.line = element_line(color = "black", size = 0.5),
      # Increase the thickness of the tick marks.
      axis.ticks = element_line(color = "black", size = 0.5),
      # Increase the font size and add bold for tick labels.
      axis.text = element_text(color = "black", size = 12),
      # Increase the font size and add bold for axis titles.
      axis.title = element_text(color = "black", size = 14)
    ) +
    theme(legend.position = legend_position) +
    # Create a single combined legend
    guides(
      # Use color for the main legend
      color = guide_legend(
        override.aes = list(
          # Explicitly define each shape for each strategy to ensure correct mapping
          shape = c(
            "Random" = NA,        # No shape for Random
            "Payoff" = 16,        # Circle for Payoff
            "Proximal" = 17,      # Triangle for Proximal
            "Prestige" = 18,      # Diamond for Prestige
            "Conformity" = 15     # Square for Conformity
          ),
          size = 3,               # Increase shape size in legend
          linetype = lty_map,     # Ensure line types match
          linewidth = 1
        ),
        title = "strategy",
        keywidth = unit(1.25, "cm") # Increase width to show dashed pattern
      ),
      # Hide the separate shape and linetype legends
      shape = "none",
      linetype = "none"
    )
  
  if (show_ci) {
    p <- p + geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                           width = 0.1, alpha = 0.7)
  }
  
  # Apply the appropriate y-axis scale based on log_scale_y parameter
  if (log_scale_y) {
    # For log scale, ensure all values are positive
    if (min(agg_data$DV_val, na.rm = TRUE) <= 0) {
      warning("Log scale requested but data contains zero or negative values. Adding small constant to make all values positive.")
      # Find the smallest positive value and use a fraction of it as offset
      min_positive <- min(agg_data$DV_val[agg_data$DV_val > 0], na.rm = TRUE)
      offset <- min_positive / 10
      # Add offset to all values
      agg_data$DV_val <- agg_data$DV_val + offset
      # Also adjust confidence intervals
      if (show_ci) {
        agg_data$ci_lower <- agg_data$ci_lower + offset
        agg_data$ci_upper <- agg_data$ci_upper + offset
      }
      # Update the plot data
      p$data <- agg_data
    }
    
    # Apply log scale with human-friendly breaks and labels showing actual numbers
    # Get the y-axis range
    y_min <- min(agg_data$DV_val, na.rm = TRUE)
    y_max <- max(agg_data$DV_val, na.rm = TRUE)
    
    # Create more comprehensive custom breaks that ensure lower values are included
    custom_breaks <- numeric(0)
    
    # Include 0 if there are any 0 values (or very small values close to 0)
    if (y_min < 0.1) {
      custom_breaks <- c(custom_breaks, 0)
    }
    
    # Add single digits (ensure we have enough coverage for lower range)
    if (y_max > 1) {
      single_digits <- c(1, 2, 3, 5)
      custom_breaks <- c(custom_breaks, single_digits[single_digits >= max(1, y_min) & single_digits <= min(9, y_max)])
    }
    
    # For the 10-100 range, include all multiples of 10 as we want both major and minor ticks
    if (y_max >= 10) {
      tens <- seq(10, min(100, y_max), by = 10)
      custom_breaks <- c(custom_breaks, tens)
    }
    
    # Add 100+ values
    if (y_max > 100) {
      higher_breaks <- c(100, 150, 200, 300, 500, 1000)
      custom_breaks <- c(custom_breaks, higher_breaks[higher_breaks <= y_max & higher_breaks > 100])
    }
    
    # If custom breaks is still empty (unlikely), fall back to automatic breaks
    if (length(custom_breaks) == 0) {
      custom_breaks <- scales::breaks_extended(n = 8)(c(y_min, y_max))
    }
    
    # Instead of using minor ticks (which might be overridden),
    # we'll include all multiples of 10 as major ticks but make some less prominent
    p <- p + scale_y_log10(
      breaks = custom_breaks,
      labels = scales::label_number()
    )
  } else if (DV == "step_payoff"){
    # Always set y-axis to 0-1.5 with ticks every 0.25, regardless of the data
    p <- p + scale_y_continuous(
      breaks = seq(0, 1.5, by = 0.5),  # Breaks from 0 to 1.5 by 0.5
      limits = c(0, 1.5)                # Hard limit at exactly 1.5
    )
  }
  
  if (!is.null(title)) {
    p <- p + labs(title = title) +  theme(plot.title = element_text(color = "black", size = 8, face = "bold", hjust = 0.5))
  }
  
  
  # Return the generated bins and x-positions along with the plot for reference
  attr(p, "bins") <- bins
  attr(p, "xposs") <- xposs
  
  if (show_plot) print(p)
  return(p)
}

plotDVbyIVBinnedRelative <- function(data, DV, DV_label, IV, IV_label,
                                     lambda_ratio = NULL,
                                     num_bins = 7,
                                     DV_trans = identity,
                                     log_scale_y = FALSE,
                                     y_range = NULL,
                                     DV_scale = NULL,
                                     # Optional bins and x-axis positions that will be generated dynamically
                                     bins = NULL,
                                     xposs = NULL,
                                     title = NULL,
                                     show_plot = FALSE,
                                     # Parameters for confidence intervals
                                     show_ci = FALSE,
                                     conf_level = 0.95,
                                     show_ylab = TRUE,
                                     legend_position = "none",
                                     relative_type = "ratio", # Options: "ratio" or "difference"
                                     # Option to show the Random baseline or not
                                     show_random = FALSE) {
  
  # Hardcoded strategy order and associated slope values.
  # (Random = 0, all others = 2)
  strategyOrder <- c("Random", "Payoff", "Proximal", "Prestige", "Conformity")
  selectedSlopes <- c(0, 2, 2, 2, 2)
  # Optionally filter by lambda_ratio. This determines which time step to choose
  # so that the time step is proportional to the size of the tree.
  if (!is.null(lambda_ratio)) {
    data <- data %>% filter(steps == floor(lambda_ratio * num_nodes))
  }
  
  if (!is.null(DV_scale)) {
    data[[DV]] <- data[[DV]] / (DV_scale * data$num_nodes)
  }
  
  if (is.null(bins)) {
    min_iv <- min(data[[IV]], na.rm = TRUE)
    max_iv <- max(data[[IV]], na.rm = TRUE)
    
    # Small epsilon for narrow min/max bins
    epsilon <- 0.001
    
    # Create equidistant x-positions
    x_positions <- seq(min_iv, max_iv, length.out = num_bins)
    
    # Create bin boundaries with narrow min/max bins
    bins <- c(min_iv - epsilon, min_iv + epsilon)
    
    # Add middle bin boundaries
    if (num_bins > 2) {
      for (i in 2:(num_bins-1)) {
        middle_point <- (x_positions[i] + x_positions[i+1]) / 2
        bins <- c(bins, middle_point)
      }
    }
    
    # Add narrow max bin
    bins <- c(bins, max_iv - epsilon, max_iv + epsilon)
    
    # Set x-positions
    xposs <- x_positions
  }
  # Initialize a data frame to store aggregated results.
  agg_data <- data.frame()
  
  # Loop over each strategy.
  for (i in seq_along(strategyOrder)) {
    strat <- strategyOrder[i]
    slope_val <- selectedSlopes[i]
    
    # Subset the data for the current strategy and slope.
    strat_data <- data %>% 
      filter(strategy == strat, slope == slope_val)
    
    # Process each bin.
    for (j in seq_len(length(bins) - 1)) {
      bin_lower <- bins[j]
      bin_upper <- bins[j + 1]
      
      # Subset rows where the IV falls inside the current bin interval.
      bin_data <- strat_data %>% 
        filter(.data[[IV]] > bin_lower, .data[[IV]] <= bin_upper)
      
      if (nrow(bin_data) == 0) next
      
      # For each unique network in the bin, choose the DV value 
      # from the middle observation.
      networks <- unique(bin_data$adj_mat)
      mid_vals <- c()
      for (net in networks) {
        net_data <- bin_data %>% filter(adj_mat == net)
        if (nrow(net_data) == 0) next
        mid_index <- ceiling(nrow(net_data) / 2)
        mid_vals <- c(mid_vals, DV_trans(net_data[[DV]][mid_index]))
      }
      
      if (length(mid_vals) > 0) {
        avg_val <- mean(mid_vals, na.rm = TRUE)
        
        # Calculate confidence intervals
        if (length(mid_vals) >= 2) {
          # Use t-distribution for confidence interval
          t_val <- qt((1 + conf_level) / 2, df = length(mid_vals) - 1)
          sd_val <- sd(mid_vals, na.rm = TRUE)
          se_val <- sd_val / sqrt(length(mid_vals))
          ci_lower <- avg_val - t_val * se_val
          ci_upper <- avg_val + t_val * se_val
        } else {
          # If only one observation, can't calculate CI
          ci_lower <- NA
          ci_upper <- NA
        }
        
        # Use the provided x position for this bin.
        agg_data <- rbind(agg_data, data.frame(
          strategy = strat,
          bin = j,
          x = xposs[j],
          DV_val = avg_val,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          n = length(mid_vals)
        ))
      }
    }
  }
  
  # For legend ordering, we want the legend to show:
  # "Payoff", "Proximal", "Prestige", "Conformity", "Random"
  agg_data$strategy <- factor(agg_data$strategy,
                              levels = c("Payoff", "Proximal", "Prestige", "Conformity", "Random"))
  
  # Create a version of the data with values relative to Random
  rel_data <- data.frame()
  
  # Process each bin to calculate relative performance
  for (j in unique(agg_data$bin)) {
    # Get all strategies for this bin
    bin_data <- agg_data %>% filter(bin == j)
    
    # Get the Random strategy value for this bin
    random_val <- bin_data %>% 
      filter(strategy == "Random") %>% 
      pull(DV_val)
    
    # Only proceed if we have a Random value for this bin
    if (length(random_val) > 0) {
      # For each strategy, calculate relative performance
      for (strat in unique(bin_data$strategy)) {
        strat_data <- bin_data %>% filter(strategy == strat)
        
        if (strat == "Random" && !show_random) {
          # Skip Random if not showing it
          next
        }
        
        # Calculate relative value
        if (relative_type == "ratio") {
          rel_val <- strat_data$DV_val / random_val
          # Also adjust confidence intervals if present
          if (show_ci) {
            rel_ci_lower <- strat_data$ci_lower / random_val
            rel_ci_upper <- strat_data$ci_upper / random_val
          }
        } else if (relative_type == "difference") {
          rel_val <- strat_data$DV_val - random_val
          # Also adjust confidence intervals if present
          if (show_ci) {
            rel_ci_lower <- strat_data$ci_lower - random_val
            rel_ci_upper <- strat_data$ci_upper - random_val
          }
        }
        
        # Add to relative data
        if (show_ci) {
          rel_data <- rbind(rel_data, data.frame(
            strategy = strat,
            bin = j,
            x = strat_data$x,
            DV_val = rel_val,
            ci_lower = rel_ci_lower,
            ci_upper = rel_ci_upper,
            n = strat_data$n
          ))
        } else {
          rel_data <- rbind(rel_data, data.frame(
            strategy = strat,
            bin = j,
            x = strat_data$x,
            DV_val = rel_val,
            n = strat_data$n
          ))
        }
      }
    }
  }
  
  # If showing Random, make sure it has the correct relative value
  if (show_random) {
    random_data <- rel_data %>% filter(strategy == "Random")
    if (relative_type == "ratio") {
      random_data$DV_val <- 1.0  # Ratio of 1.0
      if (show_ci) {
        random_data$ci_lower <- 1.0
        random_data$ci_upper <- 1.0
      }
    } else if (relative_type == "difference") {
      random_data$DV_val <- 0.0  # Difference of 0.0
      if (show_ci) {
        random_data$ci_lower <- 0.0
        random_data$ci_upper <- 0.0
      }
    }
    # Update Random in the data
    rel_data <- rel_data %>% filter(strategy != "Random")
    rel_data <- rbind(rel_data, random_data)
  }
  
  # Define manual mappings:
  col_map <- c(
    "Random" = "grey30",
    "Payoff" = "#006328",
    "Proximal" = "#ff8954",
    "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
    "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
  )
  shape_map <- c(
    "Payoff" = 16,
    "Proximal" = 17,
    "Prestige" = 18,
    "Conformity" = 15,
    "Random" = NA  # No shape for Random
  )
  lty_map <- c(
    "Random" = 2,
    "Payoff" = 1,
    "Proximal" = 1,
    "Prestige" = 1,
    "Conformity" = 1
  )
  
  shape_size_map <- c(
    "Random" = 3,
    "Payoff" = 3,  
    "Proximal" = 3, 
    "Prestige" = 3,
    "Conformity" = 3
  )
  
  
  # Create the ggplot.
  p <- ggplot(rel_data, aes(x = x, y = DV_val,
                            group = strategy,
                            color = strategy,
                            shape = strategy,
                            linetype = strategy)) +
    # Use the size mapping for line thickness
    geom_line(size = 1) 
  
  # Create a subset of data for non-Random strategies (for points)
  non_random_data <- subset(rel_data, strategy != "Random")
  
  # Add white backing points first (larger size for proper coverage)
  p <- p + geom_point(
    data = non_random_data,
    aes(size = strategy),
    color = "white",
    fill = "white",
    stroke = 2  # White border/stroke for better coverage
  )
  
  # Then add the actual colored points on top
  p <- p + geom_point(
    data = non_random_data,
    aes(size = strategy, color = strategy, shape = strategy)
  )
  
  # Apply themes and scales
  p <- p + theme_classic() +
    scale_color_manual(
      values = col_map[unique(rel_data$strategy)],
      breaks = unique(rel_data$strategy)
    ) +
    scale_shape_manual(
      values = shape_map[unique(rel_data$strategy)],
      breaks = setdiff(unique(rel_data$strategy), "Random")
    ) +
    scale_linetype_manual(
      values = lty_map[unique(rel_data$strategy)],
      breaks = unique(rel_data$strategy)
    ) +
    scale_size_manual(
      values = shape_size_map[unique(rel_data$strategy)],
      breaks = unique(rel_data$strategy),
      guide = "none"  # Hide the size legend since we only want one combined legend
    ) +
    labs(x = IV_label, 
         y = DV_label) +
    scale_x_continuous(
      limits = c(0, max(rel_data$x) * 1.05),  # Add some padding on the right
      expand = expansion(mult = c(0.05, 0.05))  # Add consistent padding on both sides
    ) +
    theme(
      # Increase the thickness of the x and y axis lines.
      axis.line = element_line(color = "black", size = 0.5),
      # Increase the thickness of the tick marks.
      axis.ticks = element_line(color = "black", size = 0.5),
      # Increase the font size and add bold for tick labels.
      axis.text = element_text(color = "black", size = 12),
      # Increase the font size and add bold for axis titles.
      axis.title = element_text(color = "black", size = 14)
    ) +
    theme(legend.position = legend_position) +
    # Create a single combined legend
    guides(
      # Use color for the main legend
      color = guide_legend(
        override.aes = list(
          shape = shape_map[unique(rel_data$strategy)]  # Only use shapes for present strategies
        ),
        title = "strategy"
      ),
      # Hide the separate shape and linetype legends
      shape = "none",
      linetype = "none"
    )
  
  # Add a reference line for relative performance
  if (relative_type == "ratio") {
    p <- p + geom_hline(yintercept = 1, linetype = "dotted", color = "black")
  } else {
    p <- p + geom_hline(yintercept = 0, linetype = "dotted", color = "black")
  }
  
  if (show_ci && "ci_lower" %in% names(rel_data) && "ci_upper" %in% names(rel_data)) {
    p <- p + geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                           width = 0.1, alpha = 0.7)
  }
  
  # Apply the appropriate y-axis scale based on log_scale_y parameter and y_range
  if (!is.null(y_range)) {
    # If y_range is provided, use it to manually set the Y-axis range
    if (log_scale_y && relative_type == "ratio") {
      p <- p + scale_y_log10(
        limits = y_range,
        breaks = scales::breaks_log(n = 6)
      )
    } else {
      p <- p + scale_y_continuous(
        limits = y_range,
        breaks = scales::pretty_breaks(n = 6)
      )
    }
  } else if (log_scale_y && relative_type == "ratio") {
    # For log scale, ensure all values are positive
    if (min(rel_data$DV_val, na.rm = TRUE) <= 0) {
      warning("Log scale requested but data contains zero or negative values. Adding small constant to make all values positive.")
      # Find the smallest positive value and use a fraction of it as offset
      min_positive <- min(rel_data$DV_val[rel_data$DV_val > 0], na.rm = TRUE)
      offset <- min_positive / 10
      # Add offset to all values
      rel_data$DV_val <- rel_data$DV_val + offset
      # Also adjust confidence intervals
      if (show_ci) {
        rel_data$ci_lower <- rel_data$ci_lower + offset
        rel_data$ci_upper <- rel_data$ci_upper + offset
      }
      # Update the plot data
      p$data <- rel_data
    }
    
    # Apply log scale with human-friendly breaks and labels showing actual numbers
    # Get the y-axis range
    y_min <- min(rel_data$DV_val, na.rm = TRUE)
    y_max <- max(rel_data$DV_val, na.rm = TRUE)
    
    # Create more comprehensive custom breaks that ensure lower values are included
    custom_breaks <- numeric(0)
    
    # Include 0 if there are any 0 values (or very small values close to 0)
    if (y_min < 0.1) {
      custom_breaks <- c(custom_breaks, 0)
    }
    
    # Add single digits (ensure we have enough coverage for lower range)
    if (y_max > 1) {
      single_digits <- c(0.5, 1, 2, 3, 5)
      custom_breaks <- c(custom_breaks, single_digits[single_digits >= max(0.5, y_min) & single_digits <= min(9, y_max)])
    }
    
    # For the 10-100 range, include all multiples of 10 as we want both major and minor ticks
    if (y_max >= 10) {
      tens <- seq(10, min(100, y_max), by = 10)
      custom_breaks <- c(custom_breaks, tens)
    }
    
    # Add 100+ values
    if (y_max > 100) {
      higher_breaks <- c(100, 150, 200, 300, 500, 1000)
      custom_breaks <- c(custom_breaks, higher_breaks[higher_breaks <= y_max & higher_breaks > 100])
    }
    
    # If custom breaks is still empty (unlikely), fall back to automatic breaks
    if (length(custom_breaks) == 0) {
      custom_breaks <- scales::breaks_extended(n = 8)(c(y_min, y_max))
    }
    
    # Instead of using minor ticks (which might be overridden),
    # we'll include all multiples of 10 as major ticks but make some less prominent
    p <- p + scale_y_log10(
      breaks = custom_breaks,
      labels = scales::label_number()
    )
  } else if (relative_type == "ratio") {
    # For ratio, we want to emphasize the 1.0 line (equal to Random)
    y_min <- min(rel_data$DV_val, na.rm = TRUE)
    y_max <- max(rel_data$DV_val, na.rm = TRUE)
    
    # Ensure 1.0 is included in the range
    y_min <- min(y_min, 0.95)
    y_max <- max(y_max, 1.05)
    
    # Set y-axis range with some padding
    p <- p + scale_y_continuous(
      limits = c(y_min * 0.95, y_max * 1.05),
      breaks = scales::pretty_breaks(n = 6)
    )
  } else if (relative_type == "difference") {
    # For difference, we want to emphasize the 0 line (equal to Random)
    y_min <- min(rel_data$DV_val, na.rm = TRUE)
    y_max <- max(rel_data$DV_val, na.rm = TRUE)
    
    # Ensure 0 is included in the range
    y_min <- min(y_min, -0.05)
    y_max <- max(y_max, 0.05)
    
    # Set y-axis range with some padding
    p <- p + scale_y_continuous(
      limits = c(y_min * 1.05, y_max * 1.05),  # More padding below
      breaks = scales::pretty_breaks(n = 6)
    )
  }
  
  if (!is.null(title)) {
    p <- p + labs(title = title) +  theme(plot.title = element_text(color = "black", size = 8, face = "bold", hjust = 0.5))
  }
  
  # Return the generated bins and x-positions along with the plot for reference
  attr(p, "bins") <- bins
  attr(p, "xposs") <- xposs
  
  if (show_plot) print(p)
  return(p)
}

plotDVbyIVSlopes <- function(data, DV, DV_label, IV, IV_label,
                             strategy, slopes = c(1.25, 5),  # Default to just weak (1.25) and strong (5) bias
                             lambda_value = NULL,
                             DV_trans = identity,
                             DV_scale = NULL,
                             bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001),
                             xposs = 2:8/2,
                             y_bins = 0.25,
                             auto_y_scale = FALSE,
                             show_title = FALSE,
                             show_ci = TRUE,
                             conf_level = 0.95,
                             legend_position = "none"
) {
  
  # Optionally filter by lambda_value.
  if (!is.null(lambda_value)) {
    data <- data %>% filter(steps == lambda_value)
  }
  
  if (!is.null(DV_scale)) {
    data[[DV]] <- data[[DV]] / DV_scale
  }
  
  # Initialize a data frame to store aggregated results.
  agg_data <- data.frame()
  
  # Get color for the strategy
  col_map <- c(
    "Random" = "grey30",
    "Payoff" = "#006328",
    "Proximal" = "#ff8954",
    "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
    "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
  )
  
  # Sort slopes to ensure correct shape assignment
  slopes <- sort(slopes)
  
  # Process data for the selected strategy and each slope value
  for (i in seq_along(slopes)) {
    slope_val <- slopes[i]
    
    # Subset the data for the current strategy and slope.
    strat_data <- data %>% 
      filter(strategy == !!strategy, slope == slope_val)
    
    # Process each bin.
    for (j in seq_len(length(bins) - 1)) {
      bin_lower <- bins[j]
      bin_upper <- bins[j + 1]
      
      # Subset rows where the IV falls inside the current bin interval.
      bin_data <- strat_data %>%
        filter(.data[[IV]] > bin_lower, .data[[IV]] <= bin_upper)
      
      if (nrow(bin_data) == 0) next
      
      # For each unique network in the bin, choose the DV value 
      # from the middle observation.
      networks <- unique(bin_data$adj_mat)
      mid_vals <- c()
      for (net in networks) {
        net_data <- bin_data %>% filter(adj_mat == net)
        if (nrow(net_data) == 0) next
        mid_index <- ceiling(nrow(net_data) / 2)
        mid_vals <- c(mid_vals, DV_trans(net_data[[DV]][mid_index]))
      }
      
      if (length(mid_vals) > 0) {
        avg_val <- mean(mid_vals, na.rm = TRUE)
        
        # Calculate confidence intervals
        if (length(mid_vals) >= 2) {
          # Use t-distribution for confidence interval
          t_val <- qt((1 + conf_level) / 2, df = length(mid_vals) - 1)
          sd_val <- sd(mid_vals, na.rm = TRUE)
          se_val <- sd_val / sqrt(length(mid_vals))
          ci_lower <- avg_val - t_val * se_val
          ci_upper <- avg_val + t_val * se_val
        } else {
          ci_lower <- NA
          ci_upper <- NA
        }
        
        # Use the provided x position for this bin.
        bias_label <- if (slope_val == 1.25) {
          "Weak bias"
        } else if (slope_val == 5) {
          "Strong bias"
        } else {
          paste("Slope =", slope_val)
        }
        
        agg_data <- rbind(agg_data, data.frame(
          slope = as.character(slope_val),
          slope_label = bias_label,
          bin = j,
          x = xposs[j],
          DV_val = avg_val,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          n = length(mid_vals)
        ))
      }
    }
  }
  
  # Also get data for Random strategy with slope = 0 for the reference line
  random_data <- data.frame()
  
  strat_data <- data %>% filter(strategy == "Random", slope == 0)
  
  # Process each bin for Random strategy
  for (j in seq_len(length(bins) - 1)) {
    bin_lower <- bins[j]
    bin_upper <- bins[j + 1]
    
    bin_data <- strat_data %>% 
      filter(.data[[IV]] > bin_lower, .data[[IV]] <= bin_upper)
    
    if (nrow(bin_data) == 0) next
    
    networks <- unique(bin_data$adj_mat)
    mid_vals <- c()
    for (net in networks) {
      net_data <- bin_data %>% filter(adj_mat == net)
      if (nrow(net_data) == 0) next
      mid_index <- ceiling(nrow(net_data) / 2)
      mid_vals <- c(mid_vals, DV_trans(net_data[[DV]][mid_index]))
    }
    
    if (length(mid_vals) > 0) {
      avg_val <- mean(mid_vals, na.rm = TRUE)
      random_data <- rbind(random_data, data.frame(
        x = xposs[j],
        DV_val = avg_val
      ))
    }
  }
  
  # Ensure we have data to plot
  if (nrow(agg_data) == 0) {
    stop("No data available for the specified strategy and slopes")
  }
  
  # Convert slope_label to factor with desired order.
  agg_data$slope_label <- factor(agg_data$slope_label, 
                                 levels = c("Weak bias", "Strong bias"))
  
  # Create the main plot
  p <- ggplot() +
    # Add Random strategy reference line if data exists
    {if (nrow(random_data) > 0)
      geom_line(
        data = random_data,
        aes(x = x, y = DV_val),
        color = "grey30",
        linetype = "dashed",
        size = 1
      )
    } +
    # Draw lines for the selected strategy (using strategy color)
    geom_line(
      data = agg_data,
      aes(
        x = x,
        y = DV_val,
        group = slope_label
      ),
      color = col_map[strategy],
      size = 1
    ) +
    # Confidence interval error bars
    {if (show_ci)
      geom_errorbar(
        data = agg_data,
        aes(
          x = x,
          ymin = ci_lower,
          ymax = ci_upper,
          group = slope_label
        ),
        color = col_map[strategy],
        width = 0.1,
        alpha = 0.7
      )
    } +
    # Draw points with shapes determined by bias label.
    geom_point(
      data = agg_data,
      aes(x = x, y = DV_val, shape = slope_label),
      color = col_map[strategy],
      size = 3
    ) +
    theme_classic() +
    scale_shape_manual(
      values = c("Weak bias" = 15, "Strong bias" = 17),
      name = "Bias strength"
    ) +
    labs(x = IV_label, y = DV_label) +
    theme(
      axis.line = element_line(color = "black", size = 1.0),
      axis.ticks = element_line(color = "black", size = 1.0),
      axis.text = element_text(color = "black", size = 12, face = "bold"),
      axis.title = element_text(color = "black", size = 14, face = "bold"),
      plot.title = element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "none"  # Hide the original legend
    ) +
    scale_y_continuous(
      breaks = seq(0, 1.5, by = 0.5),
      limits = if (!auto_y_scale) c(0, 1.5) else NULL
    ) +
    if (show_title) labs(title = strategy) else labs(title = NULL)
  
  # Create a separate data frame for the legend
  legend_data <- data.frame(
    x = 1,
    y = 1,
    slope_label = factor(c("Weak bias", "Strong bias"), 
                         levels = c("Weak bias", "Strong bias"))
  )
  
  # Create a separate plot for the legend with black shapes
  legend_plot <- ggplot(legend_data, aes(x = x, y = y, shape = slope_label)) +
    geom_point(color = "black", size = 3) +
    scale_shape_manual(
      values = c("Weak bias" = 15, "Strong bias" = 17),
      name = "Bias strength"
    ) +
    theme_void() +
    theme(legend.position = legend_position)
  
  # Extract the legend from the legend plot
  legend <- cowplot::get_legend(legend_plot)
  

  
  
  
  if (legend_position != "none") {
    grid::grid.newpage()
    grid::grid.draw(legend)
    return(legend)
  }
  print(p)
  return(p)
}


plotDVbyIVSlopesRelative <- function(data, DV, DV_label, IV, IV_label,
                                     strategy, 
                                     slopes = c(1.25, 5),  # Default to just weak (1.25) and strong (5) bias
                                     lambda_value = NULL,
                                     DV_trans = identity,
                                     DV_scale = NULL,
                                     bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001),
                                     xposs = 2:8/2,
                                     y_bins = 0.25,
                                     y_range = NULL,
                                     auto_y_scale = FALSE,
                                     show_title = FALSE,
                                     show_ci = TRUE,
                                     conf_level = 0.95,
                                     legend_position = "none",
                                     relative_type = "ratio", # Options: "ratio" or "difference"
                                     show_random = FALSE) {
  
  # Optionally filter by lambda_value.
  if (!is.null(lambda_value)) {
    data <- data %>% filter(steps == lambda_value)
  }
  
  if (!is.null(DV_scale)) {
    data[[DV]] <- data[[DV]] / DV_scale
  }
  
  # Initialize a data frame to store aggregated results.
  agg_data <- data.frame()
  
  # Get color for the strategy
  col_map <- c(
    "Random" = "grey30",
    "Payoff" = "#006328",
    "Proximal" = "#ff8954",
    "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
    "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
  )
  
  # Define shape map based on strategy, no longer based on weak/strong bias
  shape_map <- c(
    "Payoff" = 16,       # Circle
    "Proximal" = 17,     # Triangle
    "Prestige" = 18,     # Diamond
    "Conformity" = 15    # Square
  )
  
  # Define corresponding hollow shapes
  hollow_shape_map <- c(
    "Payoff" = 21,       # Hollow circle
    "Proximal" = 24,     # Hollow triangle
    "Prestige" = 23,     # Hollow diamond
    "Conformity" = 22    # Hollow square
  )
  
  # Sort slopes to ensure correct processing
  slopes <- sort(slopes)
  
  # Process data for the selected strategy and each slope value
  for (i in seq_along(slopes)) {
    slope_val <- slopes[i]
    
    # Subset the data for the current strategy and slope.
    strat_data <- data %>% 
      filter(strategy == !!strategy, slope == slope_val)
    
    # Process each bin.
    for (j in seq_len(length(bins) - 1)) {
      bin_lower <- bins[j]
      bin_upper <- bins[j + 1]
      
      # Subset rows where the IV falls inside the current bin interval.
      bin_data <- strat_data %>%
        filter(.data[[IV]] > bin_lower, .data[[IV]] <= bin_upper)
      
      if (nrow(bin_data) == 0) next
      
      # For each unique network in the bin, choose the DV value 
      # from the middle observation.
      networks <- unique(bin_data$adj_mat)
      mid_vals <- c()
      for (net in networks) {
        net_data <- bin_data %>% filter(adj_mat == net)
        if (nrow(net_data) == 0) next
        mid_index <- ceiling(nrow(net_data) / 2)
        mid_vals <- c(mid_vals, DV_trans(net_data[[DV]][mid_index]))
      }
      
      if (length(mid_vals) > 0) {
        avg_val <- mean(mid_vals, na.rm = TRUE)
        
        # Calculate confidence intervals
        if (length(mid_vals) >= 2) {
          # Use t-distribution for confidence interval
          t_val <- qt((1 + conf_level) / 2, df = length(mid_vals) - 1)
          sd_val <- sd(mid_vals, na.rm = TRUE)
          se_val <- sd_val / sqrt(length(mid_vals))
          ci_lower <- avg_val - t_val * se_val
          ci_upper <- avg_val + t_val * se_val
        } else {
          ci_lower <- NA
          ci_upper <- NA
        }
        
        # Use the provided x position for this bin.
        bias_label <- if (slope_val == 1.25) {
          "Weak bias"
        } else if (slope_val == 5) {
          "Strong bias"
        } else {
          paste("Slope =", slope_val)
        }
        
        agg_data <- rbind(agg_data, data.frame(
          slope = as.character(slope_val),
          slope_label = bias_label,
          bin = j,
          x = xposs[j],
          DV_val = avg_val,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          n = length(mid_vals)
        ))
      }
    }
  }
  
  # Get data for Random strategy with slope = 0 for the reference
  random_data <- data.frame()
  
  strat_data <- data %>% filter(strategy == "Random", slope == 0)
  
  # Process each bin for Random strategy
  for (j in seq_len(length(bins) - 1)) {
    bin_lower <- bins[j]
    bin_upper <- bins[j + 1]
    
    bin_data <- strat_data %>% 
      filter(.data[[IV]] > bin_lower, .data[[IV]] <= bin_upper)
    
    if (nrow(bin_data) == 0) next
    
    networks <- unique(bin_data$adj_mat)
    mid_vals <- c()
    for (net in networks) {
      net_data <- bin_data %>% filter(adj_mat == net)
      if (nrow(net_data) == 0) next
      mid_index <- ceiling(nrow(net_data) / 2)
      mid_vals <- c(mid_vals, DV_trans(net_data[[DV]][mid_index]))
    }
    
    if (length(mid_vals) > 0) {
      avg_val <- mean(mid_vals, na.rm = TRUE)
      random_data <- rbind(random_data, data.frame(
        bin = j,
        x = xposs[j],
        DV_val = avg_val
      ))
    }
  }
  
  # Ensure we have data to plot
  if (nrow(agg_data) == 0) {
    stop("No data available for the specified strategy and slopes")
  }
  
  # Convert slope_label to factor with desired order.
  agg_data$slope_label <- factor(agg_data$slope_label, 
                                 levels = c("Weak bias", "Strong bias"))
  
  # Create a version of the data with values relative to Random
  rel_data <- data.frame()
  
  # Process each bin to calculate relative performance
  for (j in unique(agg_data$bin)) {
    # Get the Random strategy value for this bin
    random_val <- random_data %>% 
      filter(bin == j) %>% 
      pull(DV_val)
    
    # Only proceed if we have a Random value for this bin
    if (length(random_val) > 0) {
      # For each slope in this bin, calculate relative performance
      bin_data <- agg_data %>% filter(bin == j)
      
      for (slope_label in unique(bin_data$slope_label)) {
        slope_data <- bin_data %>% filter(slope_label == slope_label)
        
        # Calculate relative value
        if (relative_type == "ratio") {
          rel_val <- slope_data$DV_val / random_val
          # Also adjust confidence intervals if present
          if (show_ci) {
            rel_ci_lower <- slope_data$ci_lower / random_val
            rel_ci_upper <- slope_data$ci_upper / random_val
          }
        } else if (relative_type == "difference") {
          rel_val <- slope_data$DV_val - random_val
          # Also adjust confidence intervals if present
          if (show_ci) {
            rel_ci_lower <- slope_data$ci_lower - random_val
            rel_ci_upper <- slope_data$ci_upper - random_val
          }
        }
        
        # Add to relative data
        if (show_ci) {
          rel_data <- rbind(rel_data, data.frame(
            slope = slope_data$slope,
            slope_label = slope_data$slope_label,
            bin = j,
            x = slope_data$x,
            DV_val = rel_val,
            ci_lower = rel_ci_lower,
            ci_upper = rel_ci_upper,
            n = slope_data$n
          ))
        } else {
          rel_data <- rbind(rel_data, data.frame(
            slope = slope_data$slope,
            slope_label = slope_data$slope_label,
            bin = j,
            x = slope_data$x,
            DV_val = rel_val,
            n = slope_data$n
          ))
        }
      }
    }
  }
  
  # Add Random baseline if needed
  if (show_random) {
    random_baseline <- random_data
    if (relative_type == "ratio") {
      random_baseline$DV_val <- 1.0  # Ratio of 1.0
    } else {
      random_baseline$DV_val <- 0.0  # Difference of 0.0
    }
    random_baseline$slope_label <- "Random baseline"
  }
  
  # Create the main plot
  p <- ggplot() +
    # Add Random strategy reference line if showing
    {if (show_random)
      geom_line(
        data = random_baseline,
        aes(x = x, y = DV_val),
        color = "grey30",
        linetype = "dashed",
        size = 1
      )
    } +
    # Draw lines for the selected strategy (using strategy color)
    geom_line(
      data = rel_data,
      aes(
        x = x,
        y = DV_val,
        group = slope_label
      ),
      color = col_map[strategy],
      size = 1
    ) +
    # Confidence interval error bars
    {if (show_ci)
      geom_errorbar(
        data = rel_data,
        aes(
          x = x,
          ymin = ci_lower,
          ymax = ci_upper,
          group = slope_label
        ),
        color = col_map[strategy],
        width = 0.1,
        alpha = 0.7
      )
    }
  
  # First add white backing shapes to hide lines
  p <- p + 
    geom_point(
      data = rel_data,
      aes(x = x, y = DV_val),
      shape = 16,  # Circle shape for backing
      color = "white",
      size = 4     # Slightly larger than the points to fully cover lines
    )
  
  # Add shape points based on strategy and bias strength
  # Use the shape determined by strategy
  # Strong bias = filled shape, Weak bias = hollow shape (outlined)
  for (bias_type in c("Strong bias", "Weak bias")) {
    bias_data <- rel_data %>% filter(slope_label == bias_type)
    
    if (nrow(bias_data) > 0) {
      if (bias_type == "Strong bias") {
        # Filled shape for strong bias
        p <- p + geom_point(
          data = bias_data,
          aes(x = x, y = DV_val),
          shape = shape_map[strategy],  # Use strategy shape (filled)
          color = col_map[strategy],
          fill = col_map[strategy],  # Fill with the same color
          size = 3
        )
      } else {
        # Hollow shape for weak bias - use filled=NA for transparent fill
        p <- p + geom_point(
          data = bias_data,
          aes(x = x, y = DV_val),
          shape = hollow_shape_map[strategy],  # Use correct hollow shape for strategy
          color = col_map[strategy],  # Border color
          fill = NA,  # Transparent fill
          size = 3,
          stroke = 1  # Border thickness
        )
      }
    }
  }
  
  # Add theme and styling
  p <- p +
    theme_classic() +
    labs(x = IV_label, y = DV_label) +
    theme(
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 14),
      plot.title = element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "none"  # Hide the original legend
    )
  
  # Add a reference line for relative performance
  if (relative_type == "ratio") {
    p <- p + geom_hline(yintercept = 1, linetype = "dotted", color = "black")
  } else {
    p <- p + geom_hline(yintercept = 0, linetype = "dotted", color = "black")
  }
  
  # Scale adjustments based on relative type and auto_y_scale
  if (!auto_y_scale) {
    if (relative_type == "ratio") {
      # For ratio, set reasonable limits that emphasize the 1.0 baseline
      p <- p + scale_y_continuous(
        breaks = seq(0, 2.5, by = 0.5),
        limits = c(0, 2.8)
      )
    } else {
      # For difference, set limits that emphasize the 0 baseline
      p <- p + scale_y_continuous(
        breaks = seq(-1, 1, by = 0.5),
        limits = c(-1, 1)
      )
    }
  } else {
    # For auto scaling, still ensure reference value is included
    y_min <- min(rel_data$DV_val, na.rm = TRUE)
    y_max <- max(rel_data$DV_val, na.rm = TRUE)
    
    if (relative_type == "ratio") {
      # Ensure 1.0 is in the range for ratio
      y_min <- min(y_min, 0.95)
      y_max <- max(y_max, 1.05)
      p <- p + scale_y_continuous(
        limits = ifelse(is.null(y_range), c(y_min * 0.95, y_max * 1.05), y_range),
        breaks = scales::pretty_breaks(n = 6)
      )
    } else {
      # Ensure 0 is in the range for difference
      y_min <- min(y_min, -0.05)
      y_max <- max(y_max, 0.05)
      p <- p + scale_y_continuous(
        limits = c(y_min * 1.05, y_max * 1.05),
        breaks = scales::pretty_breaks(n = 6)
      )
    }
  }
  
  # Add title if requested
  if (show_title) {
    if (relative_type == "ratio") {
      p <- p + labs(title = paste(strategy, "/ Random"))
    } else {
      p <- p + labs(title = paste(strategy, "- Random"))
    }
  }
  
  # Create a separate data frame for the legend
  legend_data <- data.frame(
    x = rep(1, 2),
    y = rep(1, 2),
    bias_label = factor(c("Weak bias", "Strong bias"), 
                        levels = c("Weak bias", "Strong bias"))
  )
  
  # Create a separate plot for the legend
  if (legend_position != "none") {
    # For weak bias - we need to create the "hollow" effect in the legend
    legend_plot <- ggplot() +
      # First draw points with shapes and colors
      geom_point(
        data = data.frame(x = 1, y = 1),
        aes(x = x, y = y),
        shape = shape_map[strategy],
        color = "black",
        size = 3
      ) +
      # Then for weak bias, add white center
      geom_point(
        data = data.frame(x = 1, y = 2),
        aes(x = x, y = y),
        shape = shape_map[strategy],
        color = "black",
        size = 3
      ) +
      geom_point(
        data = data.frame(x = 1, y = 2),
        aes(x = x, y = y),
        shape = 16,
        color = "white",
        size = 2
      ) +
      # Add legend labels
      scale_y_continuous(
        breaks = c(1, 2),
        labels = c("Strong bias", "Weak bias")
      ) +
      theme_void() +
      theme(
        legend.position = legend_position,
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )
    
    # Extract the legend from the legend plot
    legend <- cowplot::get_legend(legend_plot)
    grid::grid.newpage()
    grid::grid.draw(legend)
    return(legend)
  }
  
  return(p)
}



plotDVbyIV_outdeg <- function(data, DV, DV_label, IV, IV_label, lambda_value, strategy_colors = NULL) {
  if (!is.null(lambda_value)) {
    data <- data %>% filter(steps == lambda_value)
  }
  
  graph_ids <- data.frame(
    graph = unique(data$adj_mat),
    ID = 1:length(unique(data$adj_mat))
  )
  
  average_data <- data %>%
    # mutate(
    #   across(all_of(DV), ~scales::rescale(.x, to = c(0, 1)))
    # ) %>%
    group_by(adj_mat, strategy, !!sym(IV)) %>%
    summarize(avg_DV = mean(!!sym(DV), na.rm = TRUE), .groups = 'drop')
  
  if (length(unique(data$num_nodes)) > 1) {
    num_nodes <- paste0(min(data$num_nodes), " - ", max(data$num_nodes))
  } else {
    num_nodes <- data$num_nodes[1]
  }
  
  average_data$graph_id <- graph_ids$ID[match(average_data$adj_mat, graph_ids$graph)]
  average_data <- add_graph_measure(average_data, ecount, "root_outdegree")
  average_data$root_outdegree <- as.factor(average_data$root_outdegree)
  plot <- ggplot(average_data, aes_string(x = IV, y = "avg_DV", color = "root_outdegree")) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      x = IV_label,
      y = DV_label
    ) +
    theme_minimal() + geom_text(aes(label = graph_id), hjust = -0.2)
  
  # if (!is.null(strategy_colors)) {
  #   plot <- plot + scale_color_manual(name = "Strategy", values = strategy_colors)
  # } else {
  #   plot <- plot + scale_color_discrete(name = "Strategy")
  # }
  
  print(plot)
  return(plot)
}

plotDVbyIVnofilter <- function(data, DV, DV_label, IV, IV_label, strategy_colors = NULL) {
  average_data <- data %>%
    mutate(
      !!sym(DV) := scales::rescale(!!sym(DV), to = c(0, 1))
      #!!sym(IV) := scales::rescale(!!sym(IV), to = c(0, 1))
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


plotDVbyTime <- function(data, DV, DV_label, IV, IV_label, strategy_colors = NULL) {
  average_data <- data %>%
    group_by(adj_mat, strategy, !!sym(IV), structure) %>%
    summarize(avg_DV = mean(!!sym(DV), na.rm = TRUE), .groups = 'drop')
  
  plot <- ggplot(average_data, aes_string(x = IV, y = "avg_DV", color = "strategy")) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      title = paste0(DV_label, " by ", IV_label),
      x = IV_label,
      y = DV_label
    ) +
    theme_minimal() +
    facet_wrap(~ structure) 
  
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
  V(graph)$color <- rep("black", num_nodes)
  plot(graph,
       layout = layout_as_tree(graph),
       vertex.label = NA,
       vertex.size = 5,
       edge.arrow.size = 0.3
       )
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

plot_graph_panel <- function(df, output_file = "output.pdf") {
  unique_adj_matrices <- df %>% 
    distinct(adj_mat)
  
  plot_list <- list()
  
  max_plots <- ifelse(nrow(unique_adj_matrices) > 50, 50, nrow(unique_adj_matrices))
  
  for (i in 1:max_plots) {
    adj_string <- unique_adj_matrices$adj_mat[i]
    num_nodes <- sqrt(nchar(adj_string))

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
  
  combined_plot <- plot_grid(plotlist = plot_list, ncol = n_cols, nrow = n_rows, align = 'none')
  
  # Calculate PDF dimensions in inches (8K resolution: 7680x4320 pixels at 300 DPI)
  pdf_width <- 7680 / 300  # Width in inches
  pdf_height <- 4320 / 300 # Height in inches
  
  # Save the combined plot as a high-resolution PDF
  ggsave(output_file, combined_plot, width = pdf_width, height = pdf_height, dpi = 300)
  
  print(combined_plot)
  combined_plot
}


plot_variance_decomposition <- function(model_data) {
  strategies <- unique(model_data$strategy)
  results <- data.frame()
  
  for (strat in strategies) {
    subset_data <- model_data[model_data$strategy == strat, ]
    
    # Full model
    full_model <- lm(step_payoff ~ mean_prereq + root_outdeg, data = subset_data)
    r2_full <- summary(full_model)$r.squared
    
    # Calculate unique contributions
    unique_mean_prereq <- r2_full - summary(lm(step_payoff ~ root_outdeg, data = subset_data))$r.squared
    unique_root_outdeg <- r2_full - summary(lm(step_payoff ~ mean_prereq, data = subset_data))$r.squared
    shared <- r2_full - unique_mean_prereq - unique_root_outdeg
    unexplained <- 1 - r2_full
    
    results <- rbind(results, data.frame(
      strategy = strat,
      root_outdeg = unique_root_outdeg,
      mean_prereq = unique_mean_prereq,
      shared = shared,
      unexplained = unexplained
    ))
  }
  
  # Order strategies by total explained variance
  results <- results %>%
    mutate(total_explained = root_outdeg + mean_prereq + shared) %>%
    arrange(desc(total_explained))
  
  results$strategy <- factor(results$strategy, levels = results$strategy)
  
  # Reshape data for proper legend
  plot_data <- tidyr::pivot_longer(
    results,
    cols = c("root_outdeg", "mean_prereq", "shared", "unexplained"),
    names_to = "component",
    values_to = "value"
  )
  
  # Set factor levels for stacking order
  plot_data$component <- factor(plot_data$component, 
                                levels = c("root_outdeg", "mean_prereq", "shared", "unexplained"))
  
  # Create plot
  p <- ggplot(plot_data, aes(x = strategy, y = value, fill = component)) +
    geom_col(position = "stack") +
    scale_fill_manual(
      values = c("root_outdeg" = "#4e79a7", "mean_prereq" = "#f28e2b", 
                 "shared" = "#59a14f", "unexplained" = "#e15759"),
      labels = c("Root Outdegree", "Mean Prerequisites", 
                 "Shared Variance", "Unexplained"),
      name = "Component"
    ) +
    labs(x = "Strategy", y = "Proportion of Variance") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    theme_minimal()
  
  return(list(plot = p, data = results))
}


makeLegend <- function(strategies) {
  # Set up an empty plot with no visible elements
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
       xlim = c(0, 1), ylim = c(0, 1))
  
  # Define color and shape mappings
  col_map <- c(
    "Payoff" = "#006328",
    "Proximal" = "#ff8954",
    "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
    "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
  )
  
  shape_map <- c(
    "Payoff" = 16,
    "Proximal" = 17,
    "Prestige" = 18,
    "Conformity" = 15
  )
  
  # Filter mappings to only include strategies that were passed in
  colors <- col_map[strategies]
  shapes <- shape_map[strategies]
  
  # Create the legend with lines and points
  legend("center", legend = strategies, 
         col = colors, 
         pch = shapes,
         lty = 1,  # Add a line
         bty = "n",  # No box around the legend
         cex = 1.2,  # Size of the text
         pt.cex = 1.5,  # Size of the points
         lwd = 2)  # Line width
}


plotBarsbyStrategy <- function(
    data, 
    DV,
    DV_label,
    lambda_ratio = NULL
) {
  
  if (!is.null(lambda_ratio)) {
    data <- data %>% filter(steps == floor(lambda_ratio * num_nodes))
  }
  
  col_map <- c(
    "Random" = "grey30",
    "Payoff" = "#006328",
    "Proximal" = "#ff8954",
    "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
    "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
  )
  
  # Calculate means and standard errors by strategy
  summary_data <- data %>%
    group_by(strategy) %>%
    summarize(
      mean_value = mean(!!sym(DV), na.rm = TRUE),
      se = sd(!!sym(DV), na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Check if we have data
  if (nrow(summary_data) == 0) {
    stop("No data available after filtering.")
  }
  
  # Ensure strategies are in the specific order and exist in the data
  strategies <- names(col_map)[names(col_map) %in% summary_data$strategy]
  
  if (length(strategies) == 0) {
    stop("None of the specified strategies found in the data.")
  }
  
  summary_data <- summary_data[match(strategies, summary_data$strategy), ]
  summary_data <- summary_data[!is.na(summary_data$strategy), ]
  
  # Extract data for plotting
  heights <- summary_data$mean_value
  errors <- summary_data$se
  colors <- col_map[strategies]
  
  # Ensure we have finite values for plotting
  if (all(is.na(heights))) {
    stop("No valid data to plot. All values are NA.")
  }
  
  # Replace any NAs with 0 for plotting
  heights[is.na(heights)] <- 0
  errors[is.na(errors)] <- 0
  
  # Set up the plot area with pretty y-axis
  y_min <- 0
  y_max <- max(heights + errors, na.rm = TRUE) * 1.1  # Add 10% for spacing
  
  if (!is.finite(y_max) || y_max <= 0) {
    y_max <- 1  # Default if we can't determine a proper maximum
  }
  
  y_ticks <- pretty(c(y_min, y_max), n = 6)
  
  # Define x-axis range with spacing from y-axis
  x_min <- 0.5  # Start bars from 0.5 instead of 0
  x_max <- length(strategies) + 0.5
  
  # Start the plot
  par(mar = c(2, 4, 3, 2))  # Adjust margins (bottom, left, top, right)
  plot(1:length(strategies), heights, 
       type = "n",  # No plotting yet
       xlim = c(x_min, x_max),
       ylim = c(min(y_ticks), max(y_ticks)),
       axes = FALSE,
       xlab = "",
       ylab = DV_label,
       main = "",
       frame.plot = FALSE)  # No frame around the plot
  
  # Add y-axis with pretty ticks
  axis(2, at = y_ticks, las = 1)
  
  # Draw the bars
  barwidth <- 0.7
  for (i in 1:length(strategies)) {
    rect(i - barwidth/2, 0, i + barwidth/2, heights[i], 
         col = colors[i], border = NA)  # border = NA removes outlines
  }
  
  # Add error bars
  for (i in 1:length(strategies)) {
    if (errors[i] > 0) {
      arrows(i, heights[i] - errors[i], i, heights[i] + errors[i], 
             length = 0.05, angle = 90, code = 3)
    }
  }
  
  # Add a legend in top right
  legend("topright", 
         legend = strategies, 
         fill = colors,
         border = NA,
         bty = "n",  # No box around legend
         inset = c(0.02, 0.02))  # Small inset from the margins
}


plotBarsbyStrategy <- function(
    data1, 
    data2 = NULL,
    label1 = "Dataset 1",
    label2 = "Dataset 2",
    DV,
    DV_label,
    lambda_ratio = NULL
) {
  
  # Process first dataset
  if (!is.null(lambda_ratio)) {
    data1 <- data1 %>% filter(steps == floor(lambda_ratio * num_nodes))
    if (!is.null(data2)) {
      data2 <- data2 %>% filter(steps == floor(lambda_ratio * num_nodes))
    }
  }
  
  col_map <- c(
    "Random" = "grey30",
    "Payoff" = "#006328",
    "Proximal" = "#ff8954",
    "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
    "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
  )
  
  # Calculate means and standard errors for first dataset
  summary_data1 <- data1 %>%
    group_by(strategy) %>%
    summarize(
      mean_value = mean(!!sym(DV), na.rm = TRUE),
      se = sd(!!sym(DV), na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Calculate means and standard errors for second dataset if provided
  if (!is.null(data2)) {
    summary_data2 <- data2 %>%
      group_by(strategy) %>%
      summarize(
        mean_value = mean(!!sym(DV), na.rm = TRUE),
        se = sd(!!sym(DV), na.rm = TRUE) / sqrt(n()),
        .groups = 'drop'
      )
  } else {
    summary_data2 <- NULL
  }
  
  # Check if we have data
  if (nrow(summary_data1) == 0) {
    stop("No data available in the first dataset after filtering.")
  }
  
  # Ensure strategies are in the specific order and exist in the data
  strategies <- names(col_map)[names(col_map) %in% summary_data1$strategy]
  
  if (length(strategies) == 0) {
    stop("None of the specified strategies found in the data.")
  }
  
  # Align the first dataset with our strategy order
  summary_data1 <- summary_data1[match(strategies, summary_data1$strategy), ]
  summary_data1 <- summary_data1[!is.na(summary_data1$strategy), ]
  
  # Set up data for the first dataset
  heights1 <- summary_data1$mean_value
  errors1 <- summary_data1$se
  colors1 <- col_map[strategies]
  
  # Set up data for the second dataset if provided
  if (!is.null(summary_data2)) {
    # Align the second dataset with the same strategy order
    summary_data2 <- summary_data2[match(strategies, summary_data2$strategy), ]
    summary_data2 <- summary_data2[!is.na(summary_data2$strategy), ]
    
    heights2 <- summary_data2$mean_value
    errors2 <- summary_data2$se
    # We'll use the same colors but with different styling
  }
  
  # Ensure we have finite values for plotting
  if (all(is.na(heights1))) {
    stop("No valid data to plot in the first dataset. All values are NA.")
  }
  
  # Replace any NAs with 0 for plotting
  heights1[is.na(heights1)] <- 0
  errors1[is.na(errors1)] <- 0
  
  if (!is.null(summary_data2)) {
    heights2[is.na(heights2)] <- 0
    errors2[is.na(errors2)] <- 0
  }
  
  # Set up the plot area with pretty y-axis
  y_min <- 0
  
  # Find the maximum value across both datasets
  if (!is.null(summary_data2)) {
    y_max <- max(c(heights1 + errors1, heights2 + errors2), na.rm = TRUE) * 1.1
  } else {
    y_max <- max(heights1 + errors1, na.rm = TRUE) * 1.1
  }
  
  if (!is.finite(y_max) || y_max <= 0) {
    y_max <- 1  # Default if we can't determine a proper maximum
  }
  
  y_ticks <- pretty(c(y_min, y_max), n = 6)
  
  # Define x-axis positions
  n_strategies <- length(strategies)
  
  # Adjust bar width and positioning when we have two datasets
  if (!is.null(summary_data2)) {
    barwidth <- 0.3  # Narrower bars when showing two datasets
    gap <- 0.05      # Gap between paired bars
    
    # Positions for both datasets
    x_pos1 <- 1:n_strategies - gap/2 - barwidth/2
    x_pos2 <- 1:n_strategies + gap/2 + barwidth/2
    
    x_min <- 0.5
    x_max <- n_strategies + 0.5
  } else {
    barwidth <- 0.7  # Wider bars for single dataset
    x_pos1 <- 1:n_strategies
    x_min <- 0.5
    x_max <- n_strategies + 0.5
  }
  
  # Start the plot
  par(mar = c(3, 4, 3, 2))  # Adjust margins (bottom, left, top, right)
  plot(1:n_strategies, heights1, 
       type = "n",  # No plotting yet
       xlim = c(x_min, x_max),
       ylim = c(min(y_ticks), max(y_ticks)),
       axes = FALSE,
       xlab = "",
       ylab = DV_label,
       main = "",
       frame.plot = FALSE)  # No frame around the plot
  
  # Add y-axis with pretty ticks
  axis(2, at = y_ticks, las = 1)
  
  # Add x-axis with strategy labels
  axis(1, at = 1:n_strategies, labels = strategies, tick = FALSE)
  
  # Draw the bars for the first dataset
  for (i in 1:n_strategies) {
    if (!is.null(summary_data2)) {
      # If we have two datasets, position the first dataset bars to the left
      rect(x_pos1[i] - barwidth/2, 0, x_pos1[i] + barwidth/2, heights1[i], 
           col = colors1[i], border = NA)
    } else {
      # If we have only one dataset, center the bars
      rect(i - barwidth/2, 0, i + barwidth/2, heights1[i], 
           col = colors1[i], border = NA)
    }
  }
  
  # Draw the bars for the second dataset if provided
  if (!is.null(summary_data2)) {
    for (i in 1:n_strategies) {
      # Draw the second dataset bars to the right with dashed borders
      rect(x_pos2[i] - barwidth/2, 0, x_pos2[i] + barwidth/2, heights2[i], 
           col = colors1[i], 
           border = "black",      # Add black border
           lty = 2,               # Make border dashed (using numeric code for consistency)
           lwd = 1.5)             # Slightly thicker border for visibility
    }
  }
  
  # Add error bars for the first dataset
  for (i in 1:n_strategies) {
    if (errors1[i] > 0) {
      if (!is.null(summary_data2)) {
        arrows(x_pos1[i], heights1[i] - errors1[i], x_pos1[i], heights1[i] + errors1[i], 
               length = 0.05, angle = 90, code = 3)
      } else {
        arrows(i, heights1[i] - errors1[i], i, heights1[i] + errors1[i], 
               length = 0.05, angle = 90, code = 3)
      }
    }
  }
  
  # Add error bars for the second dataset if provided
  if (!is.null(summary_data2)) {
    for (i in 1:n_strategies) {
      if (errors2[i] > 0) {
        arrows(x_pos2[i], heights2[i] - errors2[i], x_pos2[i], heights2[i] + errors2[i], 
               length = 0.05, angle = 90, code = 3, 
               lty = "dashed")  # Use dashed line for error bars
      }
    }
  }
# 
#   # First, add the strategies legend
#   legend("top",
#          legend = strategies,
#          fill = colors1,
#          border = NA,
#          bty = "n",
#          horiz = TRUE,
#          inset = c(0, 0.02))  # Add inset to lower it slightly
# 
#   # Then, create a better representation for dataset legend
#   if (!is.null(summary_data2)) {
#     # Create custom graphics for the legend
#     legend("topright",
#            legend = c(label1, label2),
#            pch = 22,  # Square symbol
#            pt.bg = c("grey50", "grey50"),  # Same fill color for both
#            pt.cex = 2,  # Make squares larger
#            col = c(NA, "black"),  # No border for first, black for second
#            lty = c(0, 2),  # No line for first, dashed for second
#            lwd = c(0, 1.5),  # Line width
#            bty = "n",
#            inset = c(0.02, 0.12))  # Place below the strategy legend
#   }
}
