## First specify the packages of interest
packages <- c("BoolNet", "gt", "graphics", "grDevices")

## Now load or install & load all
package_check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
    library(x, character.only = TRUE)
  }
)

# Load the network and tools
source("shine_functions.R")
load("trap_spaces.RData")

# The fields provide structure for the simulation results
fields <- c(
  "index",
  "jps",
  "bf", 
  "os", 
  "gd", 
  "ef",
  "pef"
)


nop <- attractors$JPs == 0 & attractors$NO_b1 == 0 & 
  attractors$NO_b2 == 0 & attractors$PGI2 == 0 & attractors$Superoxide == 0
nop <- unique(attractors[nop, "attractor"])
sbf <- sum_attractors(bf, shine, attr_info)
sef <- sum_attractors(ef, shine, attr_info)
sgd <- sum_attractors(gd, shine, attr_info)

# Certain states represent ECs that contract with functional junctions 
# These may also represent normal barrier function
cbf <- attractors[attractors$MLC == 1 & attractors$JPs == 1,]
# Some attractors have some states with stable junctions and no contraction
# and others tha are different
sbf <- attractors[attractors$MLC == 0 & attractors$JPs == 1,]
# Usinig a more permissive endothelial function definition
# where vasodilator production is sufficient
eef <- attractors[attractors$PGI2 == 1 | attractors$NO_b1 == 1 | attractors$NO_b2 == 1,]

# JPs = 1, and MLC=0, without producing vasodilators
# Note that this behavior occurs without laminar flow
three <- unique(attractors[attractors$PGI2 == 0 & attractors$NO_b1 == 0 & attractors$NO_b2 == 0 & attractors$JPs == 1 & attractors$MLC == 0, "attractor"])
reps <- 10000
p3 <- perturb_attractor(
  attr_info, # getAttractors(...)
  shine, # a network
  three[1], # number of the attractor in attr_info
  5, # the maximum number of molecules affected
  reps) # repetitions
p3s <- colSums(p3$results[,2:ncol(p3$results)])
print("Behaviors reached from attractor bf no ef:")
print(p3s)

# Attractors 28 and 29 represent ECs that have functional JPs and produce
# vasodilators (NO_b1, PGI2 and no TXA2), however, they also present 
# cytoskeletal tension (MLC = 1).
cq1 <- unique(attractors[attractors$PGI2 == 1 & (attractors$NO_b1 == 1 | attractors$NO_b2 == 1) & attractors$JPs == 1 & attractors$MLC == 1, "attractor"])
cq <- cq1[cq1 < 51]
pq <- cq1[cq1 > 50 & cq1 < 61]
pq1 <- attractors[attractors$PGI2 == 1 & (attractors$NO_b1 == 1 | attractors$NO_b2 == 1) & attractors$JPs == 1 & attractors$MLC == 0, "attractor"]
pq1 <- unique(pq1[pq1 > 50])
scq <- sum_attractors(cq, shine, attr_info)
print("Characteristics of attrctors cq")
print(scq)

# perturb cuasi-quiescent (MLC = 1)attractors to see their robustnes
# and understand where noise can lead
cq_perturbation_results <- data.frame(matrix(ncol = 8, nrow = 0))
reps <- 10000
for(p in cq){
  p1 <- perturb_attractor(
    attr_info, # getAttractors(...)
    shine, # a network
    p, # number of the attractor in attr_info
    5, # the maximum number of molecules affected
    reps) # repetitions
  p1s <- colSums(p1$results[,2:ncol(p1$results)])/reps
  p1s["n"] <- p
  cq_perturbation_results <- rbind(cq_perturbation_results, p1s)
}

colnames(cq_perturbation_results) <- names(p1s)
print(cq_perturbation_results)
print(colSums(cq_perturbation_results)/nrow(cq_perturbation_results))

# Attractors where PGI2 = 1 and MLC = 1, and in half of states NO = 1 and JPs = 1

spq <- sum_attractors(pq, shine, attr_info)
print("Characteristics of pq attrctors")
print(spq)

pq_perturbation_results <- data.frame(matrix(ncol = 8, nrow = 0))
reps <- 10000
for(p in pq){
  p1 <- perturb_attractor(
    attr_info, # getAttractors(...)
    shine, # a network
    p, # number of the attractor in attr_info
    5, # the maximum number of molecules affected
    reps) # repetitions
  p1s <- colSums(p1$results[,2:ncol(p1$results)])/reps
  p1s["n"] <- p
  pq_perturbation_results <- rbind(pq_perturbation_results, p1s)
}

colnames(pq_perturbation_results) <- names(p1s)
print(pq_perturbation_results)
spqpr <- colSums(pq_perturbation_results)/nrow(pq_perturbation_results)
print(spqpr)

# Attractors where PGI2 = 1, and in half of states NO = 1 and JPs = 1  and MLC = 0

spq1 <- sum_attractors(pq1, shine, attr_info)
print("Characteristics of pq1 attrctors")
print(spq1)

pq1_perturbation_results <- data.frame(matrix(ncol = 8, nrow = 0))
reps <- 10000
for(p in pq1){
  p1 <- perturb_attractor(
    attr_info, # getAttractors(...)
    shine, # a network
    p, # number of the attractor in attr_info
    5, # the maximum number of molecules affected
    reps) # repetitions
  p1s <- colSums(p1$results[,2:ncol(p1$results)])/reps
  p1s["n"] <- p
  pq1_perturbation_results <- rbind(pq1_perturbation_results, p1s)
}

colnames(pq1_perturbation_results) <- names(p1s)
print(pq1_perturbation_results)
spq1pr <- colSums(pq1_perturbation_results)/nrow(pq1_perturbation_results)
print(spq1pr)

# PGI2 = 1 and NO = 0 and JPs = 0

ef1 <- unique(attractors[attractors$PGI2 == 1 & attractors$NO_b1 == 0 & attractors$NO_b2 == 0 & attractors$JPs == 0, "attractor"])
ef1 <- ef1[ef1 < 51]

sef1 <- sum_attractors(ef1, shine, attr_info)
print("Characteristics of ef1 attrctors")
print(sef1)

ef1_perturbation_results <- data.frame(matrix(ncol = 8, nrow = 0))
reps <- 10000
for(p in ef1){
  p1 <- perturb_attractor(
    attr_info, # getAttractors(...)
    shine, # a network
    p, # number of the attractor in attr_info
    5, # the maximum number of molecules affected
    reps) # repetitions
  p1s <- colSums(p1$results[,2:ncol(p1$results)])/reps
  p1s["n"] <- p
  ef1_perturbation_results <- rbind(ef1_perturbation_results, p1s)
}

colnames(ef1_perturbation_results) <- names(p1s)
print(ef1_perturbation_results)
sef1pr <- colSums(ef1_perturbation_results)/nrow(ef1_perturbation_results)
print(sef1pr)

# PGI2 = 1 and NO = 1 and JPs = 0
ef2pef2 <- unique(attractors[attractors$PGI2 == 1 & (attractors$NO_b1 == 1 | attractors$NO_b2 == 1) & attractors$JPs == 0 & attractors$MLC == 1, "attractor"]) 
ef2 <- ef2pef2[ef2pef2 < 51]

sef2 <- sum_attractors(ef2, shine, attr_info)
print("Characteristics of ef2 attrctors")
print(sef2)

ef2_perturbation_results <- data.frame(matrix(ncol = 8, nrow = 0))
reps <- 10000
for(p in ef2){
  p1 <- perturb_attractor(
    attr_info, # getAttractors(...)
    shine, # a network
    p, # number of the attractor in attr_info
    5, # the maximum number of molecules affected
    reps) # repetitions
  p1s <- colSums(p1$results[,2:ncol(p1$results)])/reps
  p1s["n"] <- p
  ef2_perturbation_results <- rbind(ef2_perturbation_results, p1s)
}

colnames(ef2_perturbation_results) <- names(p1s)
print(ef2_perturbation_results)
sef2pr <- colSums(ef2_perturbation_results)/nrow(ef2_perturbation_results)
print(sef2pr)

# PGI2 = 1 and NO = 1 and JPs = 0 in some states
pef <- ef2pef2[ef2pef2 > 50]

spef <- sum_attractors(pef, shine, attr_info)
print("Characteristics of pef attrctors")
print(spef)

pef_perturbation_results <- data.frame(matrix(ncol = 8, nrow = 0))
reps <- 10000
for(p in pef){
  p1 <- perturb_attractor(
    attr_info, # getAttractors(...)
    shine, # a network
    p, # number of the attractor in attr_info
    5, # the maximum number of molecules affected
    reps) # repetitions
  p1s <- colSums(p1$results[,2:ncol(p1$results)])/reps
  p1s["n"] <- p
  pef_perturbation_results <- rbind(pef_perturbation_results, p1s)
}

colnames(pef_perturbation_results) <- names(p1s)
print(pef_perturbation_results)
spefpr <- colSums(pef_perturbation_results)/nrow(pef_perturbation_results)
print(spefpr)

# PGI2 = 0 and NO = 0 and JPs = 0 and Superoxide = 0

snop <- sum_attractors(nop, shine, attr_info)
print("Characteristics of nop attrctors")
print(snop)

nop_perturbation_results <- data.frame(matrix(ncol = 8, nrow = 0))
reps <- 10000
for(p in nop){
  p1 <- perturb_attractor(
    attr_info, # getAttractors(...)
    shine, # a network
    p, # number of the attractor in attr_info
    5, # the maximum number of molecules affected
    reps) # repetitions
  p1s <- colSums(p1$results[,2:ncol(p1$results)])/reps
  p1s["n"] <- p
  nop_perturbation_results <- rbind(nop_perturbation_results, p1s)
}

colnames(nop_perturbation_results) <- names(p1s)
print(nop_perturbation_results)
snoppr <- colSums(nop_perturbation_results)/nrow(nop_perturbation_results)
print(snoppr)
