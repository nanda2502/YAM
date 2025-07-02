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
datal <- read_all(0)


data4 <- read_all(4)%>% get_default()

data <- readRDS("data.rds")

####### Figure 1 #######
data_1 <- get_default(data) 

plotDVbyIVBinned(
  data = subset(data, strategy != "Perfect"),
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
  DV_label = "Relative Performance",
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
  show_plot = T
)

###### Figure 2 e-h ######
# Varying edge weights
data_weighted <- readRDS("data8_weighted.rds")

data0 <- subset(data_weighted, mean_prereq == 0)

p2e <- plotDVbyIVBinnedRelative(
  data = data0,
  DV = "step_payoff",
  DV_label = "Relative Performance",
  IV = "edge_weight",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  y_range = c(0, 2.25)
)

data1 <- subset(data_weighted, mean_prereq == 1)

p2f <- plotDVbyIVBinnedRelative(
  data = data1,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "edge_weight",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  show_plot = T,
  y_range = c(0, 2.25)
)

data2 <- subset(data_weighted, mean_prereq == 2)

p2g <- plotDVbyIVBinnedRelative(
  data = data2,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "edge_weight",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  y_range = c(0, 2.25)
)

data3 <- subset(data_weighted, mean_prereq == 3)

p2h <- plotDVbyIVBinnedRelative(
  data = data3,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "edge_weight",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  y_range = c(0, 2.25)
)

###### Figure 2 i-l ######
# Varying slopes
data_2il <- data %>% 
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Learnability",
         edge_weight == 1,
         lambda == 0
  )

p2i <- plotDVbyIVSlopesRelative(
  data = data_2il,
  strategy = "Payoff",
  DV = "step_payoff",
  DV_label = "Relative Performance",
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

p2j <- plotDVbyIVSlopesRelative(
  data = data_2il,
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

p2k <- plotDVbyIVSlopesRelative(
  data = data_2il,
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

p2l <- plotDVbyIVSlopesRelative(
  data = data_2il,
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


###### Figure 2 m #######
# Advanced traits are more likely to be expressed

data_2m <- get_default_slopes(data) %>% 
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Depth",
         edge_weight == 1,
         lambda == 0
  )

p2m <- plotDVbyIVBinnedRelative(
  data = data_2m,
  DV = "step_payoff",
  DV_label = "Relative Performance",
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_range = c(0, 3)
)

###### Figure 2 n #######
# High payoff traits are more likely to be expressed

data_2n <- get_default_slopes(data) %>% 
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Payoffs",
         edge_weight == 1,
         lambda == 0
  )

p2n <- plotDVbyIVBinnedRelative(
  data = data_2n,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_range = c(0, 3)
)

###### Figure 2 o #######
# Advanced traits have higher payoffs
data_2o <- get_default_slopes(data) %>% 
  filter(alpha == 1,
         payoffdist == 0,
         distribution == "Learnability",
         edge_weight == 1,
         lambda == 0
  )

p2o <- plotDVbyIVBinnedRelative(
  data = data_2o,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_range = c(0, 3)
)

###### Figure 2 p #######
# Varying the skewness of payoffs

data_2p <- get_default_slopes(data) %>% 
  filter(alpha == 0,
         payoffdist == 1,
         distribution == "Learnability",
         edge_weight == 1,
         lambda == 0
  )

p2p <- plotDVbyIVBinnedRelative(
  data = data_2p,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_range = c(0, 3)
)


##### Figure 2 q-v ######

data_2q <- get_default_slopes(data) %>%
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Learnability",
         edge_weight == 1,
         lambda == 0.5,
         strategy != "Perfect"
  )

p2q <- plotDVbyIVBinned(
  data = data_2q,
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  title = "Lambda = 0.5"
)


data_2r <- get_default_slopes(data) %>%
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Learnability",
         edge_weight == 1,
         lambda == 1.0,
         strategy != "Perfect"
  )

p2r <- plotDVbyIVBinned(
  data = data_2r,
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  title = "Lambda = 1.0"
)

data_2s <- get_default_slopes(data) %>%
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Learnability",
         edge_weight == 1,
         lambda == 1.5,
         strategy != "Perfect"
  )

p2s <- plotDVbyIVBinned(
  data = data_2s,
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  title = "Lambda = 1.5"
)

data_2t <- get_default_slopes(data) %>%
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Learnability",
         edge_weight == 1,
         lambda == 2.0,
         strategy != "Perfect"
  )

p2t <- plotDVbyIVBinned(
  data = data_2t,
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  title = "Lambda = 2.0"
)

data_2u <- get_default_slopes(data) %>%
    filter(alpha == 0,
           payoffdist == 0,
           distribution == "Learnability",
           edge_weight == 1,
           lambda == 10.0,
           strategy != "Perfect"
    )

p2u <- plotDVbyIVBinned(
    data = data_2u,
    DV = "step_payoff",
    DV_label = "Performance",
    IV = "mean_prereq",
    IV_label = "Constraints",
    lambda_ratio = (5/8),
    DV_scale = (5/8),
    DV_trans = identity,
    bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
    xposs = (2:8/2) - 1,
    title = "Lambda = 10.0"
)

data_2v <- get_default_slopes(data) %>%
    filter(alpha == 0,
           payoffdist == 0,
           distribution == "Learnability",
           edge_weight == 1,
           lambda == 100.0,
           strategy != "Perfect"
    )

p2v <- plotDVbyIVBinned(
    data = data_2v,
    DV = "step_payoff",
    DV_label = "Performance",
    IV = "mean_prereq",
    IV_label = "Constraints",
    lambda_ratio = (5/8),
    DV_scale = (5/8),
    DV_trans = identity,
    bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
    xposs = (2:8/2) - 1,
    title = "Lambda = 100.0"
)

lambdafunc <- function(lambda) {
    (1-seq(0,1,by = 0.05))^lambda
    
}

plot_list <- vector("list", 6)
lambdas <- c(seq(0.5, 2, by = 0.5), 10, 100)

for (i in 1:6) {
    x <- seq(0, 1, by = 0.05)
    y <- lambdafunc(lambdas[i])
    
    df <- data.frame(x = x, y = y)
    
    plot_list[[i]] <- ggplot(df, aes(x = x, y = y)) +
        geom_line() +
        xlab("normalized missing prereqs") +
        ylab("perceived learnability") +
        theme_classic()
}

plot_grid(
    plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],plot_list[[5]], plot_list[[6]],
  p2q, p2r, p2s, p2t,p2u, p2v,
  ncol = 6,nrow = 2
) 


######## Panel ########

plot_grid(
  p2e, p2f, p2g, p2h, 
  ncol = 4
)


combined_plot <- plot_grid(
  p2a, p2b, p2c, p2d,
  p2i, p2j, p2k, p2l,
  p2m, p2n, p2o, p2p,
  p2e, p2f, p2g, p2h, 
  ncol = 4, nrow = 4,
  # Ensure equal scaling across plots
  align = 'v',
  rel_widths = c(1.15, 1, 1, 1),
  labels = LETTERS[1:16]
)

combined_plot+ theme(plot.margin = margin(2, 3, 2, 2, "mm"))

combined_plot

##### Figure 2 by Rows ######

data8 <- read_sim(800)

data_2a <- get_default(data8)

p2a <- plotDVbyIVBinnedRelative(
  data = data_2a,
  DV = "step_payoff",
  DV_label = "Relative Performance",
  IV = "mean_prereq",
  IV_label = "Constraints",
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
  IV_label = "Constraints",
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
  IV_label = "Constraints",
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
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  num_bins = 7,
  show_ci = F,
  y_range = c(0,8),
  show_plot = T
)

plot_grid(
  p2a, p2b, p2c, p2d,
  ncol = 4,
  align = 'v'
) + theme(plot.margin = margin(2, 3, 2, 2, "mm"))


data_2il <- data %>% 
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Learnability"
  )

p2i <- plotDVbyIVSlopesRelative(
  data = data_2il,
  strategy = "Payoff",
  DV = "step_payoff",
  DV_label = "Relative Performance",
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_value = 5,
  DV_scale = 5,
  DV_trans = identity,
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  y_bins = 0.25,
  auto_y_scale = F,
  show_ci = F
)

p2j <- plotDVbyIVSlopesRelative(
  data = data_2il,
  strategy = "Prestige",
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
  auto_y_scale = F,
  show_ci = F
)

p2k <- plotDVbyIVSlopesRelative(
  data = data_2il,
  strategy = "Conformity",
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
  auto_y_scale = F,
  show_ci = F
)

p2l <- plotDVbyIVSlopesRelative(
  data = data_2il,
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
  auto_y_scale = F,
  show_ci = F
)

p2legend <- plotDVbyIVSlopesRelative(
  data = data_2il,
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
  auto_y_scale = F,
  show_ci = F,
  legend_position = "right"
)

plot_grid(
  p2i, p2j, p2k, p2l,
  ncol = 4,
  align = 'v'
) + theme(plot.margin = margin(2, 3, 2, 2, "mm"))


data_weighted <- readRDS("data8_weighted.rds")

data0 <- subset(data_weighted, mean_prereq == 0)

p2e <- plotDVbyIVBinnedRelative(
  data = data0,
  DV = "step_payoff",
  DV_label = "Relative Performance",
  IV = "edge_weight",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  y_range = c(0, 2.25)
)

data1 <- subset(data_weighted, mean_prereq == 1)

p2f <- plotDVbyIVBinnedRelative(
  data = data1,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "edge_weight",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  show_plot = T,
  y_range = c(0, 2.25)
)

data2 <- subset(data_weighted, mean_prereq == 2)

p2g <- plotDVbyIVBinnedRelative(
  data = data2,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "edge_weight",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  y_range = c(0, 2.25)
)

data3 <- subset(data_weighted, mean_prereq == 3)

p2h <- plotDVbyIVBinnedRelative(
  data = data3,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "edge_weight",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  y_range = c(0, 2.25)
)

plot_grid(
  p2e, p2f, p2g, p2h, 
  ncol = 4
) + theme(plot.margin = margin(2, 3, 2, 2, "mm"))

data_default <- get_default(data)

p2default <- plotDVbyIVBinnedRelative(
  data = data_default,
  DV = "step_payoff",
  DV_label = "Relative Performance",
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
  xposs = (2:8/2) - 1,
  show_ci = F,
  y_range = c(0, 3)
)
  

# Advanced traits are more likely to be expressed
data_2m <- get_default_slopes(data) %>% 
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Depth"
  )

p2m <- plotDVbyIVBinnedRelative(
  data = data_2m,
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

# High payoff traits are more likely to be expressed

data_2n <- get_default_slopes(data) %>% 
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Payoffs"
  )

p2n <- plotDVbyIVBinnedRelative(
  data = data_2n,
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

plot_grid(
  p2default, p2m, p2n,
  nrow = 1,
  align = 'v'
) + theme(plot.margin = margin(2, 3, 2, 2, "mm"))


# Advanced traits have higher payoffs
data_2o <- get_default_slopes(data) %>% 
  filter(alpha == 1,
         payoffdist == 0,
         distribution == "Learnability"
  )

p2o <- plotDVbyIVBinnedRelative(
  data = data_2o,
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

# Varying the skewness of payoffs

data_2p <- get_default_slopes(data) %>% 
  filter(alpha == 0,
         payoffdist == 1,
         distribution == "Learnability"
  )

p2p <- plotDVbyIVBinnedRelative(
  data = data_2p,
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

plot_grid(
  p2default, p2o, p2p,
  nrow = 1,
  align = 'v'
) + theme(plot.margin = margin(2, 3, 2, 2, "mm"))

##### Figure 2 Version B ######
data8 <- read_sim(800)

data_2a <- get_default(data8)


p2a <- plotDVbyIVBinnedRelative(
  data = data_2a,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "mean_prereq",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  num_bins = 7,
  show_ci = F,
  y_range = c(0,6)
)

data30 <- read_sim(30)

data_2b <- get_default_slopes(data30)

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
  y_range = c(0,6)
)


data_weighted <- readRDS("data8_weighted.rds")


data1 <- subset(data_weighted, mean_prereq == 1)

p2c <- plotDVbyIVBinnedRelative(
  data = data1,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "edge_weight",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  show_plot = T,
  y_range = c(0, 2.25)
)

data3 <- subset(data_weighted, mean_prereq == 3)

p2d <- plotDVbyIVBinnedRelative(
  data = data3,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "edge_weight",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  y_range = c(0, 2.25)
)

data_2e <- get_default_slopes(data) %>% 
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Depth"
  )

p2e <- plotDVbyIVBinnedRelative(
  data = data_2e,
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


data_2f <- get_default_slopes(data) %>% 
  filter(alpha == 0,
         payoffdist == 0,
         distribution == "Payoffs"
  )

p2f <- plotDVbyIVBinnedRelative(
  data = data_2f,
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

data_2g <- get_default_slopes(data) %>% 
  filter(alpha == 1,
         payoffdist == 0,
         distribution == "Learnability"
  )

p2g <- plotDVbyIVBinnedRelative(
  data = data_2g,
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

data_2h <- get_default_slopes(data) %>% 
  filter(alpha == 0,
         payoffdist == 1,
         distribution == "Learnability"
  )

p2h <- plotDVbyIVBinnedRelative(
  data = data_2h,
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

plot_grid(
  p2a, p2b, p2c, p2d,
  p2e, p2f, p2g, p2h,
  ncol = 4,
  align = 'v'
) + theme(plot.margin = margin(2, 3, 2, 2, "mm"))













data_weighted2 <- read_all(84)

data0 <- subset(data_weighted2, mean_prereq == 0)

p2e <- plotDVbyIVBinnedRelative(
  data = data0,
  DV = "step_payoff",
  DV_label = "Rel. Performance",
  IV = "edge_weight",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  y_range = c(0, 2.25)
)

data1 <- subset(data_weighted2, mean_prereq == 1)

p2f <- plotDVbyIVBinnedRelative(
  data = data1,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "edge_weight",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  show_plot = T,
  y_range = c(0, 2.25)
)

data2 <- subset(data_weighted2, mean_prereq == 2)

p2g <- plotDVbyIVBinnedRelative(
  data = data2,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "edge_weight",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  y_range = c(0, 2.25)
)

data3 <- subset(data_weighted2, mean_prereq == 3)

p2h <- plotDVbyIVBinnedRelative(
  data = data3,
  DV = "step_payoff",
  DV_label = NULL,
  IV = "edge_weight",
  IV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8),
  bins = c(0.039, 0.041, 0.28, 0.44, 0.60, 0.76, 0.92, 0.999, 1.001),
  xposs = c(0.04, 0.20, 0.36, 0.52, 0.68, 0.84, 1.00),
  show_ci = F,
  y_range = c(0, 2.25)
)

































##### Table S1 ##### 
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


####### Figure 4 ######


g_baobab <- readRDS("graphs/g_baobab.rds")

p4a <- plot_graph(g_baobab)

g_tuber <- readRDS("graphs/g_tuber.rds")

p4b <- plot_graph(g_tuber)

g_cook <- readRDS("graphs/g_cook.rds")

p4c <- plot_graph(g_cook)

g_mathtech <- readRDS("graphs/g_mathtech.rds")

p4d <- plot_graph(g_mathtech)

data_baobab <- read_sim(86)

p4e <- plotBarsbyStrategy(
  data = data_baobab,
  DV = "step_payoff",
  DV_label = "Relative Performance",
  lambda_ratio = (5/8),
  DV_scale = (5/8)
)

data_tuber <- read_sim(10)

p4f <- plotBarsbyStrategy(
  data = data_tuber,
  DV = "step_payoff",
  DV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8)
)

data_cook <- read_sim(26)

p4g <- plotBarsbyStrategy(
  data = data_cook,
  DV = "step_payoff",
  DV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8)
)

data_mathtech <- read_sim(31)

p4h <- plotBarsbyStrategy(
  data = data_mathtech,
  DV = "step_payoff",
  DV_label = NULL,
  lambda_ratio = (5/8),
  DV_scale = (5/8)
)
plot_grid(
  p4a, p4b, p4c, p4d,
  ncol = 4
)
plot_grid(
  p4e, p4f, p4g, p4h,
  ncol = 4,
  labels = LETTERS[5:8]
)

plot_grid(
  p4a, p4b, p4c, p4d,
  p4e, p4f, p4g, p4h,
  ncol = 4,
  labels = LETTERS[1:8]
)

ps1 <- plotDVbyTime(
  data = data_cook,
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "steps",
  IV_label = "Time",
  title = "Restaurant Cooks",
  strategy_colors = c(
    "Random" = "grey30",
    "Payoff" = "#006328",
    "Proximal" = "#ff8954",
    "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
    "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
  )
)

ps2 <- plotDVbyTime(
  data = data_mathtech,
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "steps",
  IV_label = "Time",
  title = "Mathematical Technicians",
  strategy_colors = c(
    "Random" = "grey30",
    "Payoff" = "#006328",
    "Proximal" = "#ff8954",
    "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
    "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
  )
)

ps3 <- plotDVbyTime(
  data = data_baobab,
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "steps",
  IV_label = "Time",
  title = "Baobab Climbing",
  strategy_colors = c(
    "Random" = "grey30",
    "Payoff" = "#006328",
    "Proximal" = "#ff8954",
    "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
    "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
  )
)

ps4 <- plotDVbyTime(
  data = data_tuber,
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "steps",
  IV_label = "Time",
  title = "Tuber Digging",
  strategy_colors = c(
    "Random" = "grey30",
    "Payoff" = "#006328",
    "Proximal" = "#ff8954",
    "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
    "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
  )
)

plot_grid(
  ps1, ps2, ps3, ps4,
  ncol = 4
)
##### Figure S1 #####
# Performance over time per strategy



ps1a <- plotDVbyTimeRvals(
  data = data_1,
  DV = "step_payoff",
  DV_label = "Performance",
  strategy = "Payoff", 
  rvals = c(0, 1, 2, 3),
  show_ci = FALSE,
  conf_level = 0.95
)


##### Figure S2 #####
# Figure 1 but with individual data points to show distributions

ps1 <- plotDVbyIV(
  data = data_1,
  DV = "step_payoff",
  DV_label = "Performance",
  IV = "mean_prereq",
  IV_label = "Constraints",
  lambda_ratio = (5/8),
  DV_scale = (5/8)
)

##### Cultural Variation ####

data_var <- readRDS("data_var8.rds")

data_var2 <- read_all(89)

plotDVbyIV(
  data = data_var2,
  DV = "stationary_variation",
  DV_label = "Cultural Variation",
  IV = "mean_prereq",
  IV_label = "Constraints"
)

##### Time panel 1-10#####
# Panel showing the data from Figure 1 but at time steps 1 - 10

listS1 <- vector("list", length = 10)

for (step in 1:10) {
  
  DV_label <- if(step %% 5 == 1) "Performance" else NULL
  IV_label <- if(step > 5) "Constraints" else NULL
  
  
  listS1[[step]] <- plotDVbyIVBinned(
    data = data_1,
    DV = "step_payoff",
    DV_label = DV_label,
    IV = "mean_prereq",
    IV_label = IV_label,
    lambda_ratio = (step/8),
    DV_scale = (step/8),
    log_scale_y = F,
    DV_trans = identity,
    bins = c(0.999, 1.001, 1.25 + 1:4/2, 3.999, 4.001) - 1,
    xposs = (2:8/2) - 1,
    show_ci = F,
    legend_position = "none",
    show_plot = F,
    title = paste0("t=", step)
  )
}

plot_grid(
  listS1[[1]], listS1[[2]], listS1[[3]], listS1[[4]],
  listS1[[5]], listS1[[6]], listS1[[7]], listS1[[8]],
  listS1[[9]], listS1[[10]],
  ncol = 5
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
    group_by(mean_prereq, edge_weight, strategy) %>%
    summarize(avg_payoff = mean(step_payoff), .groups = "drop")
  
  # Find the maximum payoff for each coordinate
  max_payoffs <- payoff_data %>%
    group_by(mean_prereq, edge_weight) %>%
    summarize(max_payoff = max(avg_payoff), .groups = "drop")
  
  # Join the data and find all strategies that match the maximum payoff
  tied_strategies <- payoff_data %>%
    inner_join(max_payoffs, by = c("mean_prereq", "edge_weight")) %>%
    filter(abs(avg_payoff - max_payoff) < 1e-10)
  
  # For each coordinate, select the strategy with highest priority
  best_strategies <- tied_strategies %>%
    mutate(priority = strategy_priority[strategy]) %>%
    group_by(mean_prereq, edge_weight) %>%
    slice_min(order_by = priority, n = 1) %>%
    ungroup() %>%
    dplyr::select(-priority, -max_payoff)
  
  # Create a finer grid for visualization
  prereq_vals <- sort(unique(data$mean_prereq))
  weight_vals <- sort(unique(data$edge_weight))
  
  # Create fine grid with increased resolution for smoother appearance
  grid_size_x <- 150  
  grid_size_y <- 150
  grid_x <- seq(min(prereq_vals), max(prereq_vals), length.out = grid_size_x)
  grid_y <- seq(min(weight_vals), max(weight_vals), length.out = grid_size_y)
  
  # Use expand.grid to create all combinations
  fine_grid <- expand.grid(mean_prereq = grid_x, edge_weight = grid_y)
  
  # Convert strategies to numeric
  best_strategies$strategy_num <- match(best_strategies$strategy, strategy_levels)
  
  # Add small random jitter to data points to avoid akima artifacts with duplicates
  # This is especially important for regions with many ties
  set.seed(123)  # For reproducibility
  jittered_data <- best_strategies %>%
    group_by(mean_prereq, edge_weight) %>%
    mutate(
      mean_prereq_jitter = mean_prereq + runif(n(), -0.01, 0.01) * min(diff(sort(unique(prereq_vals)))),
      edge_weight_jitter = edge_weight + runif(n(), -0.01, 0.01) * 0.05
    ) %>%
    ungroup()
  
  # Use akima for interpolation with improved smoothing parameters
  interp_result <- akima::interp(
    x = jittered_data$mean_prereq_jitter,
    y = jittered_data$edge_weight_jitter,
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
  ggplot2::ggplot(fine_grid, ggplot2::aes(x = mean_prereq, y = edge_weight, fill = strategy)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = col_map, name = "Best Strategy") +
    ggplot2::scale_x_continuous(
      breaks = prereq_vals,
      labels = function(x) sprintf("%.1f", x)
    ) +
    ggplot2::scale_y_continuous(breaks = weight_vals) +
    ggplot2::labs(
      x = "Average number of prerequisite traits",
      y = "Strength of edge weights",
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom"
    )
}

create_strategy_heatmap_gam <- function(data) {
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
  
  # Define strategy levels
  strategy_levels <- names(col_map)
  
  # Calculate average payoff for each strategy, mean_prereq, and edge_weight
  payoff_data <- data %>%
    group_by(mean_prereq, edge_weight, strategy) %>%
    summarize(avg_payoff = mean(step_payoff), .groups = "drop")
  
  # Find the maximum payoff for each coordinate
  max_payoffs <- payoff_data %>%
    group_by(mean_prereq, edge_weight) %>%
    summarize(max_payoff = max(avg_payoff), .groups = "drop")
  
  # Join the data and find all strategies that match the maximum payoff
  highest_payoff_strategies <- payoff_data %>%
    inner_join(max_payoffs, by = c("mean_prereq", "edge_weight")) %>%
    filter(abs(avg_payoff - max_payoff) < 1e-10)
  
  # Instead of selecting based on priority, randomly select one strategy per coordinate
  # This simulates what would happen without tie-breaking
  set.seed(123)  # For reproducibility
  best_strategies <- highest_payoff_strategies %>%
    group_by(mean_prereq, edge_weight) %>%
    sample_n(1) %>%  # Randomly select one strategy per coordinate
    ungroup() %>%
    dplyr::select(-max_payoff)
  
  # Create a finer grid for visualization
  prereq_vals <- sort(unique(data$mean_prereq))
  weight_vals <- sort(unique(data$edge_weight))
  
  # Create fine grid for GAM prediction
  grid_size_x <- 600
  grid_size_y <- 600
  grid_x <- seq(min(prereq_vals), max(prereq_vals), length.out = grid_size_x)
  grid_y <- seq(min(weight_vals), max(weight_vals), length.out = grid_size_y)
  
  # Use expand.grid to create all combinations
  fine_grid <- expand.grid(mean_prereq = grid_x, edge_weight = grid_y)
  
  # Convert strategies to numeric for modeling
  best_strategies$strategy_num <- as.numeric(factor(best_strategies$strategy, levels = strategy_levels))
  
  # Fit a GAM model to predict strategy based on coordinates
  # We'll use mgcv package for GAM
  library(mgcv)
  
  # Use a GAM with smooth terms for both dimensions
  # k controls the smoothness (higher = more flexible)
  gam_model <- gam(
    strategy_num ~ s(mean_prereq, edge_weight, k=30),
    data = best_strategies,
    family = gaussian()
  )
  
  # Predict strategy numbers for the fine grid
  fine_grid$strategy_num <- predict(gam_model, newdata = fine_grid, type = "response")
  
  # Convert numeric back to strategy names, ensuring valid indices
  rounded_indices <- round(pmin(pmax(fine_grid$strategy_num, 1), length(strategy_levels)))
  fine_grid$strategy <- strategy_levels[rounded_indices]
  
  # Create the plot
  ggplot2::ggplot(fine_grid, ggplot2::aes(x = mean_prereq, y = edge_weight, fill = strategy)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = col_map, name = "Best Strategy") +
    ggplot2::scale_x_continuous(
      breaks = prereq_vals,
      labels = function(x) sprintf("%.1f", x)
    ) +
    ggplot2::scale_y_continuous(breaks = weight_vals) +
    ggplot2::labs(
      x = "Average number of prerequisite traits",
      y = "Strength of edge weights"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom"
    )
}


create_strategy_levelplot_gam <- function(data) {
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(scales)
  
  data <- data %>%
    filter(steps < 11)
  
  # Only keep relevant strategies
  data <- data %>%
    filter(strategy %in% c("Payoff", "Proximal", "Random"))
  
  col_map <- c(
    "Payoff" = "#006328",
    "Proximal" = "#ff8954"
  )
  
  # Calculate average payoff for each strategy, mean_prereq, and edge_weight
  payoff_data <- data %>%
    group_by(mean_prereq, edge_weight, strategy) %>%
    summarize(avg_payoff = mean(step_payoff), .groups = "drop")
  
  # For each coordinate, get the best of Payoff/Proximal and the Random payoff
  best_strategies <- payoff_data %>%
    filter(strategy %in% c("Payoff", "Proximal")) %>%
    group_by(mean_prereq, edge_weight) %>%
    slice_max(avg_payoff, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  random_payoff <- payoff_data %>%
    filter(strategy == "Random") %>%
    select(mean_prereq, edge_weight, random_payoff = avg_payoff)
  
  # Join to get the random payoff for each coordinate
  best_strategies <- best_strategies %>%
    left_join(random_payoff, by = c("mean_prereq", "edge_weight")) %>%
    mutate(diff_from_random = avg_payoff - random_payoff)
  
  # Fit a GAM to predict both the winning strategy and the difference
  # Convert strategy to numeric for modeling
  best_strategies$strategy_num <- as.numeric(factor(best_strategies$strategy, levels = c("Payoff", "Proximal")))
  
  # Fine grid for prediction
  prereq_vals <- sort(unique(data$mean_prereq))
  weight_vals <- sort(unique(data$edge_weight))
  grid_size_x <- 600
  grid_size_y <- 600
  grid_x <- seq(min(prereq_vals), max(prereq_vals), length.out = grid_size_x)
  grid_y <- seq(min(weight_vals), max(weight_vals), length.out = grid_size_y)
  fine_grid <- expand.grid(mean_prereq = grid_x, edge_weight = grid_y)
  
  # Fit GAMs
  gam_strategy <- gam(
    strategy_num ~ s(mean_prereq, edge_weight, k = 30),
    data = best_strategies,
    family = gaussian()
  )
  gam_diff <- gam(
    diff_from_random ~ s(mean_prereq, edge_weight, k = 30),
    data = best_strategies,
    family = gaussian()
  )
  
  # Predict on grid
  fine_grid$strategy_num <- predict(gam_strategy, newdata = fine_grid, type = "response")
  fine_grid$diff_from_random <- predict(gam_diff, newdata = fine_grid, type = "response")
  
  # Convert numeric back to strategy names
  fine_grid$strategy <- c("Payoff", "Proximal")[pmin(pmax(round(fine_grid$strategy_num), 1), 2)]
  
  # Normalize difference for alpha scaling
  diff_min <- min(fine_grid$diff_from_random, na.rm = TRUE)
  diff_max <- max(fine_grid$diff_from_random, na.rm = TRUE)
  fine_grid$alpha <- (fine_grid$diff_from_random - diff_min) / (diff_max - diff_min)
  fine_grid$alpha <- pmax(pmin(fine_grid$alpha, 1), 0.1) # avoid fully transparent
  
  # Map color and alpha
  fine_grid$fill_col <- mapply(
    function(strat, a) alpha(col_map[strat], a),
    fine_grid$strategy, fine_grid$alpha
  )
  
  # Plot
  ggplot(fine_grid, aes(x = mean_prereq, y = edge_weight)) +
    geom_tile(aes(fill = fill_col), color = NA) +
    scale_fill_identity(guide = "legend", 
                        breaks = col_map, 
                        labels = names(col_map), 
                        name = "Best Strategy"
    ) +
    scale_x_continuous(
      breaks = prereq_vals,
      labels = function(x) sprintf("%.1f", x)
    ) +
    scale_y_continuous(breaks = weight_vals) +
    labs(
      x = "Average number of prerequisite traits",
      y = "Strength of edge weights",
      title = "Level Plot: Winning Strategy and Margin over Random"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom"
    ) +
    geom_contour(aes(z = diff_from_random), color = "black")
}
create_static_surface_plot <- function(
    data,
    theta = 45,
    phi = 30,
    payoff_col = "#006328",
    proximal_col = "#ff8954"
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install 'dplyr'")
  if (!requireNamespace("mgcv", quietly = TRUE)) stop("Install 'mgcv'")
  if (!requireNamespace("plot3D", quietly = TRUE)) stop("Install 'plot3D'")
  
  library(dplyr)
  library(mgcv)
  library(plot3D)
  
  # Filter and keep only relevant strategies
  data <- data %>%
    filter(steps < 11) %>%
    filter(strategy %in% c("Payoff", "Proximal", "Random"))
  
  # Calculate average payoff for each strategy, mean_prereq, and edge_weight
  payoff_data <- data %>%
    group_by(mean_prereq, edge_weight, strategy) %>%
    summarize(avg_payoff = mean(step_payoff), .groups = "drop")
  
  # For each coordinate, get the best of Payoff/Proximal and the Random payoff
  best_strategies <- payoff_data %>%
    filter(strategy %in% c("Payoff", "Proximal")) %>%
    group_by(mean_prereq, edge_weight) %>%
    slice_max(avg_payoff, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  random_payoff <- payoff_data %>%
    filter(strategy == "Random") %>%
    select(mean_prereq, edge_weight, random_payoff = avg_payoff)
  
  # Join to get the random payoff for each coordinate
  best_strategies <- best_strategies %>%
    left_join(random_payoff, by = c("mean_prereq", "edge_weight")) %>%
    mutate(diff_from_random = avg_payoff - random_payoff)
  
  # Fit GAMs to smooth the difference and strategy
  gam_diff <- mgcv::gam(
    diff_from_random ~ s(mean_prereq, edge_weight, k = 30),
    data = best_strategies,
    family = gaussian()
  )
  gam_strat <- mgcv::gam(
    as.numeric(strategy == "Proximal") ~ s(mean_prereq, edge_weight, k = 30),
    data = best_strategies,
    family = binomial()
  )
  
  # Create a fine grid for prediction
  prereq_vals <- sort(unique(data$mean_prereq))
  weight_vals <- sort(unique(data$edge_weight))
  grid_size_x <- 50  # Reduced for better performance
  grid_size_y <- 50
  
  # Create the grid points
  x_seq <- seq(min(prereq_vals), max(prereq_vals), length.out = grid_size_x)
  y_seq <- seq(min(weight_vals), max(weight_vals), length.out = grid_size_y)
  
  # Create matrices for x and y coordinates (each point in the grid)
  x_mat <- matrix(rep(x_seq, each = grid_size_y), nrow = grid_size_x, ncol = grid_size_y)
  y_mat <- matrix(rep(y_seq, times = grid_size_x), nrow = grid_size_x, ncol = grid_size_y)
  
  # Create prediction grid
  fine_grid <- expand.grid(mean_prereq = x_seq, edge_weight = y_seq)
  
  # Predict difference and strategy on the grid
  fine_grid$diff_from_random <- predict(gam_diff, newdata = fine_grid, type = "response")
  fine_grid$proximal_prob <- predict(gam_strat, newdata = fine_grid, type = "response")
  fine_grid$strategy <- ifelse(fine_grid$proximal_prob > 0.5, "Proximal", "Payoff")
  
  # Reshape for surf3D
  z_mat <- matrix(fine_grid$diff_from_random, nrow = grid_size_x, ncol = grid_size_y)
  
  # Create a numeric matrix for coloring (1 = Payoff, 2 = Proximal)
  strat_mat <- matrix(
    ifelse(fine_grid$strategy == "Payoff", 1, 2),
    nrow = grid_size_x, 
    ncol = grid_size_y
  )
  
  # Use plot3D::persp3D with matrices of the same dimensions
  plot3D::persp3D(
    x = x_mat,
    y = y_mat,
    z = z_mat,
    colvar = strat_mat,  # Color by strategy
    col = c(payoff_col, proximal_col),  # Colors for strategies
    border = "black",
    lwd = 0.1,
    shade = 0.5,
    ticktype = "detailed",
    xlab = "Average number of prerequisite traits",
    ylab = "Strength of edge weights",
    zlab = "Difference from Random",
    main = "Best Strategy (Color) and Difference from Random (Height)",
    theta = theta,
    phi = phi,
    colkey = list(
      at = c(1.25, 1.75),  # Position ticks in the middle of each color
      labels = c("Payoff", "Proximal"),
      side = 4,
      length = 0.5,
      width = 0.5
    )
  )
}









data_weighted85 <- read_all(85)

create_strategy_heatmap_gam(subset(data_weighted85, strategy != "Random" & steps == 5))

create_strategy_levelplot_gam(subset(data_weighted85, steps == 5))

level_data <- data_weighted85 %>%
  filter(steps == 5) %>%
  mutate(step_payoff = step_payoff/5)

create_static_surface_plot(level_data, theta = 210, phi = 10)


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
  plot <- ggplot(data[data$steps == 5 & data$strategy == strategy,], aes(x = mean_prereq, y = step_payoff, color = as.factor(slope), group = slope)) +
    geom_smooth(method = "loess", se = FALSE) +
    geom_smooth(data = data[data$steps == 5 & data$strategy == "Random", ],
                aes(x = avg_path_length, y = step_payoff), 
                method = "loess", se = FALSE, color = "black", linetype = "dashed") + 
    labs(x = "mean distance to root", y = "performance", color = "strength of bias") +
    ggtitle(strategy) +
    theme_minimal() + 
    ylim(min(data$step_payoff), max(data$step_payoff))
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







legend_only_plot <- function(col_map) {
  library(ggplot2)
  
  # Create a data frame with one row per color
  df <- data.frame(
    x = rep(1, length(col_map)),
    y = 1:length(col_map),
    group = factor(names(col_map), levels = names(col_map))
  )
  
  # Create the plot
  p <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_point(size = 3) +  # Add points (will be hidden later)
    scale_color_manual(values = col_map) +
    theme(
      # Remove all plot elements except legend
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.background = element_blank(),
      
      # Format the legend
      legend.position = "right",
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(shape = 16))) +  # Use circles (shape 16)
    # This hides the plot but keeps the legend
    coord_cartesian(xlim = c(0, 0), ylim = c(0, 0))
  
  return(p)
}

col_map <- c(
  "Payoff" = "#006328",
  "Proximal" = "#ff8954",
  "Prestige" = adjustcolor("#cb5b85", alpha.f = 0.5),
  "Conformity" = adjustcolor("#0163c2", alpha.f = 0.5)
)

# Create and display the legend-only plot
legend_only_plot(col_map)







# First, let's analyze how performance scales with R for each strategy
# We'll focus on performance at a specific time step to compare strategies

# Extract unique R values (mean_prereq) and strategies in the data
unique_R_values <- unique(data_1$mean_prereq)
unique_strategies <- unique(data_1$strategy)

# Look at data at step 5 (arbitrary choice to represent mid-term performance)
performance_by_R_and_strategy <- data_1 %>%
  dplyr::filter(steps == 5) %>%
  dplyr::group_by(mean_prereq, strategy) %>%
  dplyr::summarize(
    mean_performance = mean(step_payoff),
    sd_performance = sd(step_payoff),
    n = n(),
    se_performance = sd_performance / sqrt(n)
  )

# Print a summary of the results
print("Performance by R value and strategy at step 5:")
print(performance_by_R_and_strategy)

# Fit models to understand the relationship between R and performance for each strategy
# We'll try different functional forms: linear, quadratic, and exponential

fit_results <- list()
for (strat in unique_strategies) {
  data_for_strategy <- performance_by_R_and_strategy %>%
    dplyr::filter(strategy == strat)
  
  # Linear model: P = a - b*R
  linear_model <- lm(mean_performance ~ mean_prereq, data = data_for_strategy)
  
  # Quadratic model: P = a - b*R^2
  data_for_strategy$R_squared <- data_for_strategy$mean_prereq^2
  quadratic_model <- lm(mean_performance ~ R_squared, data = data_for_strategy)
  
  # Exponential model: log(P) = a - b*R (equivalent to P = exp(a) * exp(-b*R))
  exponential_model <- lm(log(mean_performance) ~ mean_prereq, data = data_for_strategy)
  
  fit_results[[strat]] <- list(
    linear = summary(linear_model),
    quadratic = summary(quadratic_model),
    exponential = summary(exponential_model)
  )
}

# Print model fit results
print("Model fitting results:")
print(fit_results)

# Find crossover points between Payoff bias and Proximal learning
# First, extract coefficients for the best-fitting models
if ("Payoff" %in% unique_strategies && "Proximal" %in% unique_strategies) {
  payoff_coeffs <- fit_results[["Payoff"]]$quadratic$coefficients
  proximal_coeffs <- fit_results[["Proximal"]]$linear$coefficients
  
  # Assuming Payoff is best fit by quadratic: a_p - b_p*R^2
  # And Proximal is best fit by linear: a_x - b_x*R
  a_p <- payoff_coeffs[1, 1]  # Intercept for payoff
  b_p <- abs(payoff_coeffs[2, 1])  # Coefficient for R^2
  
  a_x <- proximal_coeffs[1, 1]  # Intercept for proximal
  b_x <- abs(proximal_coeffs[2, 1])  # Coefficient for R
  
  # Solve: a_p - b_p*R^2 = a_x - b_x*R for R
  # Rearranging: b_p*R^2 - b_x*R + (a_x - a_p) = 0
  c <- a_x - a_p
  
  # Quadratic formula: R = (-b_x ± sqrt(b_x^2 - 4*b_p*c)) / (2*b_p)
  discriminant <- b_x^2 - 4 * b_p * c
  
  if (discriminant >= 0) {
    r_crossover1 <- (-b_x + sqrt(discriminant)) / (2 * b_p)
    r_crossover2 <- (-b_x - sqrt(discriminant)) / (2 * b_p)
    
    print(paste("Potential crossover points between Payoff and Proximal at R =", 
                round(r_crossover1, 3), "or", round(r_crossover2, 3)))
    
    # Determine which solution is valid in our R range
    valid_crossovers <- c()
    if (r_crossover1 >= min(unique_R_values) && r_crossover1 <= max(unique_R_values)) {
      valid_crossovers <- c(valid_crossovers, r_crossover1)
    }
    if (r_crossover2 >= min(unique_R_values) && r_crossover2 <= max(unique_R_values)) {
      valid_crossovers <- c(valid_crossovers, r_crossover2)
    }
    
    if (length(valid_crossovers) > 0) {
      print(paste("Valid crossover points in our R range:", paste(round(valid_crossovers, 3), collapse = ", ")))
    } else {
      print("No valid crossover points in our R range.")
    }
  } else {
    print("No real crossover points exist.")
  }
}

# Visualize the relationship between R and performance for each strategy
# Create prediction data
r_seq <- seq(min(unique_R_values), max(unique_R_values), length.out = 100)
prediction_data <- expand.grid(mean_prereq = r_seq, strategy = unique_strategies)
prediction_data$R_squared <- prediction_data$mean_prereq^2

# Add predictions for each model type
predictions <- data.frame()
for (strat in unique_strategies) {
  data_for_strategy <- prediction_data %>%
    dplyr::filter(strategy == strat)
  
  # Get coefficients from the fitted models
  linear_coef <- fit_results[[strat]]$linear$coefficients
  quadratic_coef <- fit_results[[strat]]$quadratic$coefficients
  exponential_coef <- fit_results[[strat]]$exponential$coefficients
  
  # Calculate predictions
  linear_pred <- linear_coef[1, 1] + linear_coef[2, 1] * data_for_strategy$mean_prereq
  quadratic_pred <- quadratic_coef[1, 1] + quadratic_coef[2, 1] * data_for_strategy$R_squared
  exponential_pred <- exp(exponential_coef[1, 1] + exponential_coef[2, 1] * data_for_strategy$mean_prereq)
  
  # Add to predictions dataframe
  strat_predictions <- data_for_strategy %>%
    dplyr::mutate(
      linear_pred = linear_pred,
      quadratic_pred = quadratic_pred,
      exponential_pred = exponential_pred
    )
  
  predictions <- rbind(predictions, strat_predictions)
}

# Print the top rows of predictions
print("Sample of predictions:")
print(head(predictions))

# Examine performance trajectories over time for different R values
# Focus on R=0 (no constraints) and R=3 (high constraints)
time_trajectories <- data_1 %>%
  dplyr::filter(mean_prereq %in% c(0, 3)) %>%
  dplyr::group_by(mean_prereq, strategy, steps) %>%
  dplyr::summarize(
    mean_performance = mean(step_payoff),
    sd_performance = sd(step_payoff)
  )

print("Performance trajectories over time:")
print(head(time_trajectories))

# Calculate performance ratios relative to Random strategy
performance_ratios <- data_1 %>%
  dplyr::filter(steps == 5) %>%
  dplyr::select(mean_prereq, strategy, step_payoff, adj_mat) %>%
  dplyr::group_by(mean_prereq, adj_mat) %>%
  dplyr::mutate(
    random_payoff = step_payoff[strategy == "Random"],
    ratio_to_random = step_payoff / random_payoff
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(mean_prereq, strategy) %>%
  dplyr::summarize(
    mean_ratio = mean(ratio_to_random, na.rm = TRUE),
    sd_ratio = sd(ratio_to_random, na.rm = TRUE)
  )

print("Performance ratios relative to Random strategy:")
print(performance_ratios)

