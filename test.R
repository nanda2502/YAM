calculate_variance_decomposition <- function() {
    # Set up parallel processing with all cores
    plan(multisession, workers = availableCores())
    
    # Define discrete distributions
    gap_ranges <- list(close = 1:2, intermediate = 3:4, far = 5:6)
    payoff_means <- list(high = 90, medium = 75, low = 50)
    noise_range <- -7:7
    
    configs <- list(
        A = list(gaps = c("far", "close", "intermediate"), payoffs = c("high", "medium", "low")),
        B = list(gaps = c("far", "close", "intermediate"), payoffs = c("high", "low", "medium")),
        C = list(gaps = c("intermediate", "close", "far"), payoffs = c("high", "medium", "low")),
        D = list(gaps = c("intermediate", "close", "far"), payoffs = c("high", "low", "medium"))
    )
    
    strategies <- c("gap_first", "pure_payoff", "pure_proximal", "pure_control", "mixed_50_50", "uniform_third")
    
    # Create a grid of all computation chunks
    computation_grid <- expand.grid(
        config_name = names(configs),
        constraint = c("low", "high"),
        opt1_gap = 1:2,  # indices into gap ranges
        opt2_gap = 1:2,
        opt3_gap = 1:2,
        stringsAsFactors = FALSE
    )
    
    cat(sprintf("Processing %d chunks across %d cores...\n", 
                nrow(computation_grid), availableCores()))
    
    # Process each chunk in parallel
    chunk_results <- future_map(1:nrow(computation_grid), function(i) {
        row <- computation_grid[i, ]
        config_name <- row$config_name
        constraint <- row$constraint
        config <- configs[[config_name]]
        
        payoff_opt <- which(config$payoffs == "high")
        proximal_opt <- which(config$gaps == "close")
        control_opt <- setdiff(1:3, c(payoff_opt, proximal_opt))[1]
        
        # Get actual gap values
        opt1_gap <- gap_ranges[[config$gaps[1]]][row$opt1_gap]
        opt2_gap <- gap_ranges[[config$gaps[2]]][row$opt2_gap]
        opt3_gap <- gap_ranges[[config$gaps[3]]][row$opt3_gap]
        gaps <- c(opt1_gap, opt2_gap, opt3_gap)
        
        gap_prob <- 1 / (2^3)
        
        strategy_outcomes <- list()
        for (strategy in strategies) {
            strategy_outcomes[[strategy]] <- list(values = c(), weights = c(), 
                                                  config = c(), constraint = c(), gap = c())
        }
        
        # Process all noise combinations for this gap configuration
        for (received_noise1 in noise_range) {
            for (received_noise2 in noise_range) {
                for (received_noise3 in noise_range) {
                    noise_prob <- 1 / (length(noise_range)^3)
                    
                    payoffs <- c(
                        pmax(0, payoff_means[[config$payoffs[1]]] + received_noise1),
                        pmax(0, payoff_means[[config$payoffs[2]]] + received_noise2),
                        pmax(0, payoff_means[[config$payoffs[3]]] + received_noise3)
                    )
                    
                    base_prob <- 0.25 * 0.5 * gap_prob * noise_prob
                    
                    strategy_choices <- list(
                        pure_payoff = payoff_opt,
                        pure_proximal = proximal_opt,
                        pure_control = control_opt
                    )
                    
                    if (constraint == "low") {
                        strategy_choices$gap_first <- payoff_opt
                        
                        for (strat_name in c("gap_first", "pure_payoff", "pure_proximal", "pure_control")) {
                            chosen <- strategy_choices[[strat_name]]
                            strategy_outcomes[[strat_name]]$values <- c(strategy_outcomes[[strat_name]]$values, payoffs[chosen])
                            strategy_outcomes[[strat_name]]$weights <- c(strategy_outcomes[[strat_name]]$weights, base_prob)
                            strategy_outcomes[[strat_name]]$config <- c(strategy_outcomes[[strat_name]]$config, config_name)
                            strategy_outcomes[[strat_name]]$constraint <- c(strategy_outcomes[[strat_name]]$constraint, constraint)
                            strategy_outcomes[[strat_name]]$gap <- c(strategy_outcomes[[strat_name]]$gap, paste(gaps, collapse=","))
                        }
                        
                        strategy_outcomes$mixed_50_50$values <- c(strategy_outcomes$mixed_50_50$values, payoffs[payoff_opt], payoffs[proximal_opt])
                        strategy_outcomes$mixed_50_50$weights <- c(strategy_outcomes$mixed_50_50$weights, base_prob * 0.5, base_prob * 0.5)
                        strategy_outcomes$mixed_50_50$config <- c(strategy_outcomes$mixed_50_50$config, config_name, config_name)
                        strategy_outcomes$mixed_50_50$constraint <- c(strategy_outcomes$mixed_50_50$constraint, constraint, constraint)
                        strategy_outcomes$mixed_50_50$gap <- c(strategy_outcomes$mixed_50_50$gap, paste(gaps, collapse=","), paste(gaps, collapse=","))
                        
                        for (opt in 1:3) {
                            strategy_outcomes$uniform_third$values <- c(strategy_outcomes$uniform_third$values, payoffs[opt])
                            strategy_outcomes$uniform_third$weights <- c(strategy_outcomes$uniform_third$weights, base_prob / 3)
                            strategy_outcomes$uniform_third$config <- c(strategy_outcomes$uniform_third$config, config_name)
                            strategy_outcomes$uniform_third$constraint <- c(strategy_outcomes$uniform_third$constraint, constraint)
                            strategy_outcomes$uniform_third$gap <- c(strategy_outcomes$uniform_third$gap, paste(gaps, collapse=","))
                        }
                        
                    } else {
                        if (min(gaps) == 1) {
                            strategy_choices$gap_first <- which.min(gaps)
                        } else {
                            expected_payoffs <- payoffs / gaps
                            strategy_choices$gap_first <- which.max(expected_payoffs)
                        }
                        
                        for (strat_name in c("gap_first", "pure_payoff", "pure_proximal", "pure_control")) {
                            chosen <- strategy_choices[[strat_name]]
                            success_prob <- 1 / gaps[chosen]
                            
                            strategy_outcomes[[strat_name]]$values <- c(strategy_outcomes[[strat_name]]$values, payoffs[chosen], 0)
                            strategy_outcomes[[strat_name]]$weights <- c(strategy_outcomes[[strat_name]]$weights, base_prob * success_prob, base_prob * (1 - success_prob))
                            strategy_outcomes[[strat_name]]$config <- c(strategy_outcomes[[strat_name]]$config, config_name, config_name)
                            strategy_outcomes[[strat_name]]$constraint <- c(strategy_outcomes[[strat_name]]$constraint, constraint, constraint)
                            strategy_outcomes[[strat_name]]$gap <- c(strategy_outcomes[[strat_name]]$gap, paste(gaps, collapse=","), paste(gaps, collapse=","))
                        }
                        
                        for (opt in c(payoff_opt, proximal_opt)) {
                            success_prob <- 1 / gaps[opt]
                            strategy_outcomes$mixed_50_50$values <- c(strategy_outcomes$mixed_50_50$values, payoffs[opt], 0)
                            strategy_outcomes$mixed_50_50$weights <- c(strategy_outcomes$mixed_50_50$weights, base_prob * 0.5 * success_prob, base_prob * 0.5 * (1 - success_prob))
                            strategy_outcomes$mixed_50_50$config <- c(strategy_outcomes$mixed_50_50$config, config_name, config_name)
                            strategy_outcomes$mixed_50_50$constraint <- c(strategy_outcomes$mixed_50_50$constraint, constraint, constraint)
                            strategy_outcomes$mixed_50_50$gap <- c(strategy_outcomes$mixed_50_50$gap, paste(gaps, collapse=","), paste(gaps, collapse=","))
                        }
                        
                        for (opt in 1:3) {
                            success_prob <- 1 / gaps[opt]
                            strategy_outcomes$uniform_third$values <- c(strategy_outcomes$uniform_third$values, payoffs[opt], 0)
                            strategy_outcomes$uniform_third$weights <- c(strategy_outcomes$uniform_third$weights, base_prob / 3 * success_prob, base_prob / 3 * (1 - success_prob))
                            strategy_outcomes$uniform_third$config <- c(strategy_outcomes$uniform_third$config, config_name, config_name)
                            strategy_outcomes$uniform_third$constraint <- c(strategy_outcomes$uniform_third$constraint, constraint, constraint)
                            strategy_outcomes$uniform_third$gap <- c(strategy_outcomes$uniform_third$gap, paste(gaps, collapse=","), paste(gaps, collapse=","))
                        }
                    }
                }
            }
        }
        
        return(strategy_outcomes)
    }, .options = furrr_options(seed = TRUE), .progress = TRUE)
    
    # Combine results from all chunks
    combined_outcomes <- list()
    for (strategy in strategies) {
        combined_outcomes[[strategy]] <- list(values = c(), weights = c(), config = c(), constraint = c(), gap = c())
        for (chunk_result in chunk_results) {
            combined_outcomes[[strategy]]$values <- c(combined_outcomes[[strategy]]$values, chunk_result[[strategy]]$values)
            combined_outcomes[[strategy]]$weights <- c(combined_outcomes[[strategy]]$weights, chunk_result[[strategy]]$weights)
            combined_outcomes[[strategy]]$config <- c(combined_outcomes[[strategy]]$config, chunk_result[[strategy]]$config)
            combined_outcomes[[strategy]]$constraint <- c(combined_outcomes[[strategy]]$constraint, chunk_result[[strategy]]$constraint)
            combined_outcomes[[strategy]]$gap <- c(combined_outcomes[[strategy]]$gap, chunk_result[[strategy]]$gap)
        }
    }
    
    # Calculate weighted means and variances
    cat("\nStrategy Performance:\n")
    cat(sprintf("%-20s %10s %10s\n", "Strategy", "Mean", "Variance"))
    cat(strrep("-", 42), "\n")
    
    results <- list()
    for (strategy in strategies) {
        values <- combined_outcomes[[strategy]]$values
        weights <- combined_outcomes[[strategy]]$weights
        
        weighted_mean <- sum(values * weights) / sum(weights)
        weighted_var <- sum(weights * (values - weighted_mean)^2) / sum(weights)
        
        cat(sprintf("%-20s %10.2f %10.2f\n", strategy, weighted_mean, weighted_var))
        
        results[[strategy]] <- list(
            mean = weighted_mean,
            variance = weighted_var,
            outcomes = combined_outcomes[[strategy]]
        )
    }
    
    # Variance decomposition
    cat("\n\nVariance Decomposition:\n")
    cat(strrep("=", 80), "\n")
    
    for (strategy in strategies) {
        values <- results[[strategy]]$outcomes$values
        weights <- results[[strategy]]$outcomes$weights
        configs <- results[[strategy]]$outcomes$config
        constraints <- results[[strategy]]$outcomes$constraint
        gaps <- results[[strategy]]$outcomes$gap
        
        total_mean <- results[[strategy]]$mean
        total_var <- results[[strategy]]$variance
        
        config_means <- tapply(1:length(values), configs, function(idx) {
            sum(values[idx] * weights[idx]) / sum(weights[idx])
        })
        between_config_var <- sum(sapply(names(config_means), function(cfg) {
            idx <- which(configs == cfg)
            sum(weights[idx]) * (config_means[cfg] - total_mean)^2
        })) / sum(weights)
        
        constraint_means <- tapply(1:length(values), constraints, function(idx) {
            sum(values[idx] * weights[idx]) / sum(weights[idx])
        })
        between_constraint_var <- sum(sapply(names(constraint_means), function(cons) {
            idx <- which(constraints == cons)
            sum(weights[idx]) * (constraint_means[cons] - total_mean)^2
        })) / sum(weights)
        
        gap_means <- tapply(1:length(values), gaps, function(idx) {
            sum(values[idx] * weights[idx]) / sum(weights[idx])
        })
        between_gap_var <- sum(sapply(names(gap_means), function(g) {
            idx <- which(gaps == g)
            sum(weights[idx]) * (gap_means[g] - total_mean)^2
        })) / sum(weights)
        
        within_gap_var <- total_var - between_gap_var
        
        cat(sprintf("\n%s:\n", toupper(strategy)))
        cat(sprintf("  Total variance:              %10.2f\n", total_var))
        cat(sprintf("  Due to configuration:        %10.2f (%5.1f%%)\n", 
                    between_config_var, 100 * between_config_var / total_var))
        cat(sprintf("  Due to constraint condition: %10.2f (%5.1f%%)\n", 
                    between_constraint_var, 100 * between_constraint_var / total_var))
        cat(sprintf("  Due to gap realization:      %10.2f (%5.1f%%)\n", 
                    between_gap_var, 100 * between_gap_var / total_var))
        cat(sprintf("  Residual (noise + success):  %10.2f (%5.1f%%)\n", 
                    within_gap_var, 100 * within_gap_var / total_var))
    }
    
    return(results)
}

results <- calculate_variance_decomposition()