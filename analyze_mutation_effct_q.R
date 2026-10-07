# Load the network and tools
source("shine_functions.R")
shine <- loadNetwork("nets\\levels.bnet")

# Analyze the wild type
initial_states <- 1000
reps <- 10

na <- analyze_net(shine, initial_states, reps)
tssv <- na$tssv
tssv_sd <- na$tssv_sd

# Analyze mutants

a <- plot_mutant_behavior(shine, c('Gs'), c(0), 'Gs_0', tssv, tssv_sd)
b <- plot_mutant_behavior(shine, c('Gs'), c(1), 'Gs_1', tssv, tssv_sd)
rGs <- list(all, a, b)
d_Gs <- plot_node_behavior('Gs', rGs)
c <- plot_mutant_behavior(shine, c('AC'), c(0), 'AC_0', tssv, tssv_sd)
d <- plot_mutant_behavior(shine, c('AC'), c(1), 'AC_1', tssv, tssv_sd)
rAC <- list(all, c, d)
d_AC <- plot_node_behavior('AC', rAC)
e <- plot_mutant_behavior(shine, c('PGI2'), c(0), 'PGI2_0', tssv, tssv_sd)
f <- plot_mutant_behavior(shine, c('PGI2'), c(1), 'PGI2_1', tssv, tssv_sd)
rPGI2 <- list(all, e, f)
d_PGI2 <- plot_node_behavior('PGI2', rPGI2)
g <- plot_mutant_behavior(shine, c('TIE2'), c(0), 'TIE2_0', tssv, tssv_sd)
h <- plot_mutant_behavior(shine, c('TIE2'), c(1), 'TIE2_1', tssv, tssv_sd)
rTIE2 <- list(all, g, h)
d_TIE2 <- plot_node_behavior('TIE2', rTIE2)
i <- plot_mutant_behavior(shine, c('JPs'), c(0), 'JPs_0', tssv, tssv_sd)
j <- plot_mutant_behavior(shine, c('JPs'), c(1), 'JPs_1', tssv, tssv_sd)
rJPs <- list(all, i, j)
d_JPs <- plot_node_behavior('JPs', rJPs)
k <- plot_mutant_behavior(shine, c('RAC1'), c(0), 'RAC1_0', tssv, tssv_sd)
l <- plot_mutant_behavior(shine, c('RAC1'), c(1), 'RAC1_1', tssv, tssv_sd)
rRAC1 <- list(all, k, l)
d_JPs <- plot_node_behavior('RAC1', rRAC1)