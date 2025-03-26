##### Load Packages #####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(igraph)
library(combinat)
library(ggplot2)
library(gridExtra)
library(ggraph)
library(grid)
library(cowplot)
library(lattice)
library(ppcor)

source("preprocessing.R")
source("plotting.R")

##### Data ######
data <- read_all(8)
saveRDS(data, "data.rds")
data <- readRDS("data.rds")

####### Figure 1 #######
data_1 <- get_default(data)

plotDVbyIVBinned(
  data = data_1,
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  log_scale_y = F,
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  show_ci = F,
  legend_position = "right",
  show_plot = T
)
plotDVbyIVBinnedRelative(
  data = data_1,
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "mean_prereq",
  IV_label = "Average number of prerequisite traits (R)",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  log_scale_y = F,
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  show_ci = F,
  show_random = F
)


##### Figure 2 ########

###### Figure 2 a-d ######
# Larger structures

data8 <- read_sim(800)

data_2a <- get_default(data8)

p2a <- plotDVbyIVBinnedRelative(
  data = data_2a,
  DV = "step_payoff",
  DV_label = "Rel. Performance",
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  show_ci = F,
  show_plot = T,
  y_range = c(0, 8)
)

data20 <- read_sim(20)

data_2b <- get_default_slopes(data20)

p2b <- plotDVbyIVBinnedRelative(
  data = data_2b,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  num_bins = 7,
  show_ci = F,
  y_range = c(0,8),
  show_plot = T
)

data30 <- read_sim(30)

data_2c <- get_default_slopes(data30)

p2c <- plotDVbyIVBinnedRelative(
  data = data_2c,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  num_bins = 7,
  show_ci = F,
  y_range = c(0,8),
)

data50 <- read_sim(50)

data_2d <- get_default_slopes(data50)

p2d <- plotDVbyIVBinnedRelative(
  data = data_2d,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  num_bins = 7,
  show_ci = F,
  y_range = c(0,8),
  show_plot = T,
)

###### Figure 2 e-h ######
# Varying slopes
data_2eh <- data %>% 
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Learnability"
  )

p2e <- plotDVbyIVSlopesRelative(
  data = data_2eh,
  strategy = "Payoff",
  DV = "step_payoff",
  DV_label = "Rel. Performance",
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_value = 5,
  DV_scale = 5,
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_bins = 0.25,
  auto_y_scale = F,
  show_ci = F
)

p2f <- plotDVbyIVSlopesRelative(
  data = data_2eh,
  strategy = "Prestige",
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_value = 5,
  DV_scale = 5,
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_bins = 0.25,
  auto_y_scale = F,
  show_ci = F
)

p2g <- plotDVbyIVSlopesRelative(
  data = data_2eh,
  strategy = "Conformity",
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_value = 5,
  DV_scale = 5,
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_bins = 0.25,
  auto_y_scale = F,
  show_ci = F
)

p2h <- plotDVbyIVSlopesRelative(
  data = data_2eh,
  strategy = "Proximal",
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_value = 5,
  DV_scale = 5,
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_bins = 0.25,
  auto_y_scale = F,
  show_ci = F
)

p2legend <- plotDVbyIVSlopes(
  data = data_2eh,
  strategy = "Proximal",
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_value = 5,
  DV_scale = 5,
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_bins = 0.25,
  legend_position = "right"
)


###### Figure 2 i #######
# Advanced traits are more likely to be expressed

data_2i <- get_default_slopes(data) %>% 
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Depth"
  )

p2i <- plotDVbyIVBinnedRelative(
  data = data_2i,
  DV = "step_payoff",
  DV_label = "Rel. Performance",
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_range = c(0, 3)
)

###### Figure 2 j #######
# High payoff traits are more likely to be expressed

data_2j <- get_default_slopes(data) %>% 
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Payoffs"
  )

p2j <- plotDVbyIVBinnedRelative(
  data = data_2j,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_range = c(0, 3)
)

###### Figure 2 k #######
# Advanced traits have higher payoffs
data_2k <- get_default_slopes(data) %>% 
  filter(alpha == 1,
         payoffdist == 0,
         distribution == "Learnability"
  )

p2k <- plotDVbyIVBinnedRelative(
  data = data_2k,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_range = c(0, 3)
)

###### Figure 2 l #######
# Varying the skewness of payoffs

data_2l <- get_default_slopes(data) %>% 
  filter(alpha == 0,
         payoffdist == 1,
         distribution == "Learnability"
  )

p2l <- plotDVbyIVBinnedRelative(
  data = data_2l,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_range = c(0, 3)
)



######## Panel ########



combined_plot <- plot_grid(
  p2a, p2b, p2c, p2d,
  p2e, p2f, p2g, p2h, 
  p2i, p2j, p2k, p2l,
  ncol = 4, nrow = 3,
  # Ensure equal scaling across plots
  align = 'h',
  rel_widths = c(1.15, 1, 1, 1)
)

combined_plot+ theme(plot.margin = margin(0, 1, 0, 0, "mm"))

##### Fig S1 ##### 
# Variance decomposition of step_payoff


root_outdegree <- function(g) {
  degree(g, 1, mode = "out")
}

data <- add_graph_measure(data,root_outdegree, "root_outdeg")

# Prepare the dataset
model_data <- get_default(data[data$steps == 5, ])

full_model <- lm(step_payoff ~ (mean_prereq + root_outdeg) * strategy, data = model_data)

# Basic ANOVA table
anova_table <- anova(full_model)

anova_df <- as.data.frame(anova_table)

# Format the Sum Sq column specifically
anova_df$`Sum Sq` <- format(anova_df$`Sum Sq`, digits = 2, nsmall = 3)

# Print the modified data frame
print(anova_df)

# Calculate partial eta squared for each term
# Formula: SS_effect / (SS_effect + SS_residual)
SS_residual <- anova_table["Residuals", "Sum Sq"]
partial_eta_squared <- anova_table[-nrow(anova_table), "Sum Sq"] / 
  (anova_table[-nrow(anova_table), "Sum Sq"] + SS_residual)

# Add partial eta squared to the table
anova_table$`Partial eta^2` <- c(partial_eta_squared, NA)

# Print the enhanced ANOVA table
print(anova_table)

# Fit two nested models (excluding the interaction model)
# Model 1: Only mean_prereq (with strategy as control)
model1 <- lm(step_payoff ~ mean_prereq + strategy, data = model_data)

# Model 2: Both main effects (with strategy as control)
model2 <- lm(step_payoff ~ mean_prereq + root_outdeg + strategy, data = model_data)

# Calculate R-squared for each model
r2_model1 <- summary(model1)$r.squared
r2_model2 <- summary(model2)$r.squared

# Calculate the incremental R-squared
r2_increment_prereq <- r2_model1  # mean_prereq alone
r2_increment_rootdeg <- r2_model2 - r2_model1  # additional from root_outdeg

# Print R-squared values and increments
cat("R² values:\n")
cat("Model 1 (mean_prereq): ", r2_model1, "\n")
cat("Model 2 (+ root_outdeg): ", r2_model2, "\n\n")

cat("Incremental R² values:\n")
cat("mean_prereq: ", r2_increment_prereq, "\n")
cat("root_outdeg: ", r2_increment_rootdeg, "\n")

# Create data for the stacked bar plot - using actual R² values * 100 for the y-axis
# For the first bar: just mean_prereq
bar1 <- data.frame(
  Model = "Mean Prerequisites",
  Component = "Mean Prerequisites",
  Value = r2_model1 * 100  # Convert to percentage (88.4%)
)

# For the second bar: mean_prereq + root_outdeg
bar2 <- data.frame(
  Model = rep("Main Effects", 2),
  Component = c("Mean Prerequisites", "Root Outdegree"),
  Value = c(r2_model1 * 100, r2_increment_rootdeg * 100)  # Convert to percentages (88.4% + 6.0%)
)

# Combine all bars
plot_data <- rbind(bar1, bar2)

# Set factor levels for proper ordering
plot_data$Model <- factor(plot_data$Model, 
                          levels = c("Mean Prerequisites", "Main Effects"))
plot_data$Component <- factor(plot_data$Component, 
                              levels = c("Mean Prerequisites", "Root Outdegree"))

# Create the stacked bar plot
ggplot(plot_data, aes(x = Model, y = Value, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = c("Mean Prerequisites" = "#4575B4", 
                               "Root Outdegree" = "#D73027")) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(
    x = NULL,
    y = "R² Value (%)",
    fill = "Component",
    title = "Incremental Variance Explained by Model Components"
  ) +
  # Add labels showing the R² values for each model - positioned above bars
  geom_text(data = data.frame(
    Model = factor(c("Mean Prerequisites", "Main Effects"), 
                   levels = c("Mean Prerequisites", "Main Effects")),
    y = c(91.5, 97.5),
    label = c("R² = 0.884", "R² = 0.944")
  ),
  aes(x = Model, y = y, label = label),
  size = 3.5,
  inherit.aes = FALSE) +
  ylim(0, 100) +
  scale_y_continuous(breaks = seq(0, 100, by = 20))

####### Hosseinioun et al. ######
data_weighted <- read_sim(121)

data_unweighted <- read_sim(122)

plotBarsbyStrategy(
  data1 = data_weighted,
  data2 = data_unweighted,
  label1 = "Weighted",
  label2 = "Unweighted",
  DV = "step_payoff",
  DV_label = "Performance",
  lambda_ratio = 5/8
)



plotBarsbyStrategy(
  data = data_weighted,
  DV = "step_payoff",
  DV_label = "Performance",
  lambda_ratio = 5/8
)

plotDVbyIVBinned(
  data = data_real,
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "steps",
  IV_label = "Time",
  num_bins = 7,
  log_scale_y = F,
  legend_position = "none"
)


##### Heatmap #####
create_strategy_heatmap <- function(data) {
  data <- data %>%
    filter(steps < 11)
  
  # Ensure the color map is properly applied
  col_map <- c(
    "Random" = "grey30",
    "Payoff" = "#006328",
    "Proximal" = "#ff8954",
    "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
    "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
  )
  
  # Define strategy levels and priorities
  strategy_levels <- names(col_map)
  strategy_priority <- c("Payoff" = 1, "Proximal" = 2, "Prestige" = 3, "Conformity" = 4, "Random" = 5)
  
  # Calculate average payoff for each strategy, mean_prereq, and steps
  payoff_data <- data %>%
    group_by(mean_prereq, steps, strategy) %>%
    summarize(avg_payoff = mean(step_payoff), .groups = "drop")
  
  # Find the maximum payoff for each coordinate
  max_payoffs <- payoff_data %>%
    group_by(mean_prereq, steps) %>%
    summarize(max_payoff = max(avg_payoff), .groups = "drop")
  
  # Join the data and find all strategies that match the maximum payoff
  tied_strategies <- payoff_data %>%
    inner_join(max_payoffs, by = c("mean_prereq", "steps")) %>%
    filter(abs(avg_payoff - max_payoff) < 1e-10)
  
  # For each coordinate, select the strategy with highest priority
  best_strategies <- tied_strategies %>%
    mutate(priority = strategy_priority[strategy]) %>%
    group_by(mean_prereq, steps) %>%
    slice_min(order_by = priority, n = 1) %>%
    ungroup() %>%
    dplyr::select(-priority, -max_payoff)
  
  # Create a finer grid for visualization
  prereq_vals <- sort(unique(data$mean_prereq))
  steps_vals <- 1:10
  
  # Create fine grid with increased resolution for smoother appearance
  grid_size_x <- 300  # Increased from 150
  grid_size_y <- 300  # Increased from 150
  grid_x <- seq(min(prereq_vals), max(prereq_vals), length.out = grid_size_x)
  grid_y <- seq(min(steps_vals), max(steps_vals), length.out = grid_size_y)
  
  # Use expand.grid to create all combinations
  fine_grid <- expand.grid(mean_prereq = grid_x, steps = grid_y)
  
  # Convert strategies to numeric
  best_strategies$strategy_num <- match(best_strategies$strategy, strategy_levels)
  
  # Add small random jitter to data points to avoid akima artifacts with duplicates
  # This is especially important for regions with many ties
  set.seed(123)  # For reproducibility
  jittered_data <- best_strategies %>%
    group_by(mean_prereq, steps) %>%
    mutate(
      mean_prereq_jitter = mean_prereq + runif(n(), -0.01, 0.01) * min(diff(sort(unique(prereq_vals)))),
      steps_jitter = steps + runif(n(), -0.01, 0.01) * 0.05
    ) %>%
    ungroup()
  
  # Use akima for interpolation with improved smoothing parameters
  interp_result <- akima::interp(
    x = jittered_data$mean_prereq_jitter,
    y = jittered_data$steps_jitter,
    z = jittered_data$strategy_num,
    xo = grid_x,
    yo = grid_y,
    linear = TRUE,     # Try linear interpolation for smoother transitions
    extrap = TRUE,     # Allow extrapolation
    duplicate = "mean" # Handle duplicates by averaging
  )
  
  # Apply a smoothing filter to reduce jaggedness
  # Use a simple 3x3 mean filter
  smooth_z <- interp_result$z
  for (i in 2:(nrow(smooth_z)-1)) {
    for (j in 2:(ncol(smooth_z)-1)) {
      window <- smooth_z[(i-1):(i+1), (j-1):(j+1)]
      smooth_z[i, j] <- mean(window, na.rm = TRUE)
    }
  }
  interp_result$z <- smooth_z
  
  # Convert the interpolation result to a data frame
  fine_grid$strategy_num <- as.vector(interp_result$z)
  
  # Convert numeric back to strategy names, ensuring valid indices
  rounded_indices <- round(pmin(pmax(fine_grid$strategy_num, 1), length(strategy_levels)))
  fine_grid$strategy <- strategy_levels[rounded_indices]
  
  # Create the plot
  ggplot2::ggplot(fine_grid, ggplot2::aes(x = mean_prereq, y = steps, fill = strategy)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = col_map, name = "Best Strategy") +
    ggplot2::scale_x_continuous(
      breaks = prereq_vals,
      labels = function(x) sprintf("%.1f", x)
    ) +
    ggplot2::scale_y_continuous(breaks = steps_vals) +
    ggplot2::labs(
      x = "Average number of prerequisite traits",
      y = "Time constraints"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom"
    )
}

create_strategy_heatmap(subset(get_default(data), strategy != "Random"))





##### Payoff shuffles sample size #####

data <- read.csv("./output/raw_values.csv") %>% 
  filter(step == 5)

sample_sizes <- unique(c(50, 100, 200, 500, 1000, 2000, 5000))
total_shuffles <- max(data$shuffle_idx) + 1  # Add 1 because indices start at 0
sample_sizes <- sample_sizes[sample_sizes <= total_shuffles]

# Create a dataframe to store results
convergence_results <- data.frame()

# For each sample size, calculate mean and sd
for (size in sample_sizes) {
  # Take a subset of the data (first 'size' shuffles)
  subset_data <- data[data$shuffle_idx < size, ]
  
  # Calculate statistics
  mean_payoff <- mean(subset_data$payoff)
  sd_payoff <- sd(subset_data$payoff)
  cv_payoff <- (sd_payoff / mean_payoff) * 100  # Coefficient of variation in percentage
  
  # Calculate 95% confidence interval
  ci_width <- qt(0.975, df = size - 1) * sd_payoff / sqrt(size)
  ci_lower <- mean_payoff - ci_width
  ci_upper <- mean_payoff + ci_width
  relative_ci_width <- (ci_width / mean_payoff) * 100  # CI width as percentage of mean
  
  # Add to results
  result_row <- data.frame(
    sample_size = size,
    mean = mean_payoff,
    sd = sd_payoff,
    cv = cv_payoff,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    relative_ci_width = relative_ci_width
  )
  
  convergence_results <- rbind(convergence_results, result_row)
}

# Plot mean convergence by step
ggplot(convergence_results, aes(x = sample_size, y = mean)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Convergence of Mean Payoff at Step 5 with Increasing Sample Size",
       x = "Sample Size (log scale)",
       y = "Mean Payoff")


# Plot relative CI width convergence
ggplot(convergence_results, aes(x = sample_size, y = relative_ci_width)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red") +  # 5% reference line
  geom_hline(yintercept = 1, linetype = "dashed", color = "green") +  # 1% reference line
  theme_minimal() +
  labs(title = "Convergence of Relative CI Width at Step 5 with Increasing Sample Size",
       x = "Sample Size (log scale)",
       y = "Relative CI Width (%)")

# Calculate percent change in mean and SD between consecutive sample sizes
convergence_results$mean_change_pct <- c(NA, diff(convergence_results$mean) / convergence_results$mean[-nrow(convergence_results)] * 100)
convergence_results$sd_change_pct <- c(NA, diff(convergence_results$sd) / convergence_results$sd[-nrow(convergence_results)] * 100)

# Create formatted table for display
formatted_table <- convergence_results %>%
  mutate(
    mean = round(mean, 6),
    sd = round(sd, 6),
    cv = round(cv, 2),
    mean_change_pct = round(mean_change_pct, 2),
    sd_change_pct = round(sd_change_pct, 2),
    relative_ci_width = round(relative_ci_width, 2)
  ) %>%
  select(sample_size, mean, mean_change_pct, sd, sd_change_pct, cv, relative_ci_width)

# Print the table with percent changes
print(formatted_table)

# Find the smallest sample size where:
# 1. The relative CI width is < 5%
ci_threshold <- 5
ci_indices <- which(convergence_results$relative_ci_width < ci_threshold)
min_size_ci <- ifelse(length(ci_indices) > 0, 
                      convergence_results$sample_size[min(ci_indices)], 
                      Inf)

# 2. The percent change in mean is < 1%
mean_change_threshold <- 1
mean_indices <- which(abs(convergence_results$mean_change_pct) < mean_change_threshold & !is.na(convergence_results$mean_change_pct))
min_size_mean <- ifelse(length(mean_indices) > 0, 
                        convergence_results$sample_size[min(mean_indices)], 
                        Inf)

# 3. The percent change in SD is < 5%
sd_change_threshold <- 5
sd_indices <- which(abs(convergence_results$sd_change_pct) < sd_change_threshold & !is.na(convergence_results$sd_change_pct))
min_size_sd <- ifelse(length(sd_indices) > 0, 
                      convergence_results$sample_size[min(sd_indices)], 
                      Inf)

# Print the recommended sample sizes
cat("\nRecommended minimum sample sizes:\n")
cat("Based on CI width < 5%:", ifelse(is.finite(min_size_ci), min_size_ci, "Not reached"), "\n")
cat("Based on mean stability < 1%:", ifelse(is.finite(min_size_mean), min_size_mean, "Not reached"), "\n")
cat("Based on SD stability < 5%:", ifelse(is.finite(min_size_sd), min_size_sd, "Not reached"), "\n")

# Overall recommendation
recommended_size <- max(min_size_ci, min_size_mean, min_size_sd)
cat("\nOverall recommended minimum sample size:", 
    ifelse(is.finite(recommended_size), recommended_size, "Need more samples"))


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
data_perf$DV_scaled <- scales::rescale(data_perf$step_payoff, to = c(0, 1))
plot_graph_panel(data_perf[data_perf$DV_scaled > 0.6 & data_perf$DV_scaled < 0.7,])
plot_graph_panel(data_perf[data_perf$DV_scaled > 0.4 & data_perf$DV_scaled < 0.5,])

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


data <- read.csv("adj_mat_20.csv", colClasses = c(adj_mat = "character"))

plot_graph_panel(data)


data <- add_graph_measure(data, calc_avg_path_length, "avg_path_length")


hist(data$avg_path_length, breaks = 50, main = "Average path length", xlab = "Average path length")




