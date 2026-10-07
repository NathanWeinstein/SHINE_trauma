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

a <- plot_mutant_behavior(shine, c('S1P_b1', 'S1P_b2'), c(0, 0), 'S1P_0', tssv, tssv_sd)
b <- plot_mutant_behavior(shine, c('S1P_b1', 'S1P_b2'), c(1, 0), 'S1P_1', tssv, tssv_sd)
c <- plot_mutant_behavior(shine, c('S1P_b1', 'S1P_b2'), c(0, 1), 'S1P_2', tssv, tssv_sd)
rS1P <- list(all, a, b, c)
d_S1P <- plot_node_behavior('S1P', rS1P)
d <- plot_mutant_behavior(shine, c('Gi_b1', 'Gi_b2'), c(0, 0), 'Gi_0', tssv, tssv_sd)
e <- plot_mutant_behavior(shine, c('Gi_b1', 'Gi_b2'), c(1, 0), 'Gi_1', tssv, tssv_sd)
f <- plot_mutant_behavior(shine, c('Gi_b1', 'Gi_b2'), c(0, 1), 'Gi_2', tssv, tssv_sd)
rGi <- list(all, d, e, f)
d_Gi <- plot_node_behavior('Gi', rGi)
g <- plot_mutant_behavior(shine, c('NO_b1', 'NO_b2'), c(0, 0), 'NO_0', tssv, tssv_sd)
h <- plot_mutant_behavior(shine, c('NO_b1', 'NO_b2'), c(1, 0), 'NO_1', tssv, tssv_sd)
i <- plot_mutant_behavior(shine, c('NO_b1', 'NO_b2'), c(0, 1), 'NO_2', tssv, tssv_sd)
rNO <- list(all, g, h, i)
d_NO <- plot_node_behavior('NO', rNO)
j <- plot_mutant_behavior(shine, c('Ca_b1', 'Ca_b2'), c(0, 0), 'Ca_0', tssv, tssv_sd)
k <- plot_mutant_behavior(shine, c('Ca_b1', 'Ca_b2'), c(1, 0), 'Ca_1', tssv, tssv_sd)
l <- plot_mutant_behavior(shine, c('Ca_b1', 'Ca_b2'), c(0, 1), 'Ca_2', tssv, tssv_sd)
rCa <- list(all, j, k, l)
d_Ca <- plot_node_behavior('Ca', rCa)
m <- plot_mutant_behavior(shine, c('VEGFA_b1', 'VEGFA_b2'), c(0, 0), 'VEGFA_0', tssv, tssv_sd)
n <- plot_mutant_behavior(shine, c('VEGFA_b1', 'VEGFA_b2'), c(1, 0), 'VEGFA_1', tssv, tssv_sd)
o <- plot_mutant_behavior(shine, c('VEGFA_b1', 'VEGFA_b2'), c(0, 1), 'VEGFA_2', tssv, tssv_sd)
rVEGFA <- list(all, m, n, o)
d_VEGFA <- plot_node_behavior('VEGFA', rVEGFA)
