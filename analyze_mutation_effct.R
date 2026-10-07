
# Load the network and tools
source("shine_functions.R")
shine <- loadNetwork("nets\\levels.bnet")

# Analyze the wild type
initial_states <- 1000
reps <- 10

na <- analyze_net(shine, initial_states, reps)
tssv <- na$tssv
tssv_sd <- na$tssv_sd
all <- list(tssv = tssv, tssv_sd = tssv_sd)

# Analyze mutants

a <- plot_mutant_behavior(shine, c('HIF1'), c(0), 'HIF1_0', tssv, tssv_sd)
b <- plot_mutant_behavior(shine, c('HIF1'), c(1), 'HIF1_1', tssv, tssv_sd)
rhif1 <- list(all, a, b)
d_hif1 <- plot_node_behavior('HIF1', rhif1)
c <- plot_mutant_behavior(shine, c('HIF1', 'SS_b1', 'SS_b2'), c(1, 0, 0), 'HIF1_1_SS_0', tssv, tssv_sd, save = TRUE)
e <- plot_mutant_behavior(shine, c('SS_b1', 'SS_b2'), c(0, 0), 'SS_0', tssv, tssv_sd)
f <- plot_mutant_behavior(shine, c('SS_b1', 'SS_b2'), c(1, 0), 'SS_1', tssv, tssv_sd)
g <- plot_mutant_behavior(shine, c('SS_b1', 'SS_b2'), c(0, 1), 'SS_2', tssv, tssv_sd)
rSS <- list(all, e, f, g)
d_SS <- plot_node_behavior('SS', rSS)
h <- plot_mutant_behavior(shine, c('WNT7a'), c(0), 'WNT7a_0', tssv, tssv_sd)
i <- plot_mutant_behavior(shine, c('WNT7a'), c(1), 'WNT7a_1', tssv, tssv_sd)
rWNT7a <- list(all, h, i)
d_WNT7a <- plot_node_behavior('WNT7a', rWNT7a)
j <- plot_mutant_behavior(shine, c('NE_b1', 'NE_b2'), c(0, 0), 'NE_0', tssv, tssv_sd)
k <- plot_mutant_behavior(shine, c('NE_b1', 'NE_b2'), c(1, 0), 'NE_1', tssv, tssv_sd)
l <- plot_mutant_behavior(shine, c('NE_b1', 'NE_b2'), c(0, 1), 'NE_2', tssv, tssv_sd)
rNE <- list(all, j, k, l)
d_NE <- plot_node_behavior('NE', rNE)

