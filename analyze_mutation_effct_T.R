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

a <- plot_mutant_behavior(shine, c('Gq'), c(0), 'Gq_0', tssv, tssv_sd)
b <- plot_mutant_behavior(shine, c('Gq'), c(1), 'Gq_1', tssv, tssv_sd)
rGq <- list(all, a, b)
d_Gq <- plot_node_behavior('Gq', rGq)
c <- plot_mutant_behavior(shine, c('ANGII'), c(0), 'ANGII_0', tssv, tssv_sd)
d <- plot_mutant_behavior(shine, c('ANGII'), c(1), 'ANGII_1', tssv, tssv_sd)
rANGII <- list(all, c, d)
d_ANGII <- plot_node_behavior('ANGII', rANGII)
e <- plot_mutant_behavior(shine, c('G12'), c(0), 'G12_0', tssv, tssv_sd)
f <- plot_mutant_behavior(shine, c('G12'), c(1), 'G12_1', tssv, tssv_sd)
rG12 <- list(all, e, f)
d_G12 <- plot_node_behavior('G12', rG12)
g <- plot_mutant_behavior(shine, c('TXA2'), c(0), 'TXA2_0', tssv, tssv_sd)
h <- plot_mutant_behavior(shine, c('TXA2'), c(1), 'TXA2_1', tssv, tssv_sd)
rTXA2 <- list(all, g, h)
d_TXA2 <- plot_node_behavior('TXA2', rTXA2)
i <- plot_mutant_behavior(shine, c('ONOO'), c(0), 'ONOO_0', tssv, tssv_sd)
j <- plot_mutant_behavior(shine, c('ONOO'), c(1), 'ONOO_1', tssv, tssv_sd)
rONOO <- list(all, i, j)
d_ONOO <- plot_node_behavior('ONOO', rONOO)
k <- plot_mutant_behavior(shine, c('ETS'), c(0), 'ETS_0', tssv, tssv_sd)
l <- plot_mutant_behavior(shine, c('ETS'), c(1), 'ETS_1', tssv, tssv_sd)
rETS <- list(all, k, l)
d_ETS <- plot_node_behavior('ETS', rETS)
m <- plot_mutant_behavior(shine, c('Superoxide'), c(0), 'Superoxide_0', tssv, tssv_sd)
n <- plot_mutant_behavior(shine, c('Superoxide'), c(1), 'Superoxide_1', tssv, tssv_sd)
rSuperoxide <- list(all, m, n)
d_Superoxide <- plot_node_behavior('Superoxide', rSuperoxide)
o <- plot_mutant_behavior(shine, c('SRC_b1', 'SRC_b2'), c(0, 0), 'SRC_0', tssv, tssv_sd)
p <- plot_mutant_behavior(shine, c('SRC_b1', 'SRC_b2'), c(1, 0), 'SRC_1', tssv, tssv_sd)
q <- plot_mutant_behavior(shine, c('SRC_b1', 'SRC_b2'), c(0, 1), 'SRC_2', tssv, tssv_sd)
rSRC <- list(all, o, p, q)
d_SRC <- plot_node_behavior('SRC', rSRC)
r <- plot_mutant_behavior(shine, c('RHOA'), c(0), 'RHOA_0', tssv, tssv_sd, save = TRUE)
s <- plot_mutant_behavior(shine, c('RHOA'), c(1), 'RHOA_1', tssv, tssv_sd)
rRHOA <- list(all, r, s)
d_RHOA <- plot_node_behavior('RHOA', rRHOA)
t <- plot_mutant_behavior(shine, c('MLC'), c(0), 'MLC_0', tssv, tssv_sd)
u <- plot_mutant_behavior(shine, c('MLC'), c(1), 'MLC_1', tssv, tssv_sd)
rMLC <- list(all, t, u)
d_MLC <- plot_node_behavior('MLC', rMLC)

