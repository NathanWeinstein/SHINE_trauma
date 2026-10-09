## First specify the packages of interest
packages <- c("BoolNet", "deSolve")

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

if (!dir.exists("Figures\\ode")) {
  dir.create("Figures\\ode", recursive = TRUE)
}

# Load the network and tools
source("shine_functions.R")
load("trap_spaces.RData")

#  0 < hk , γk and 0 ≤ xk , ωk ≤ 1.
squad <- function(x_i, w_i, h, g_i){
  num <- -exp(0.5 * h) + exp(-h * (w_i - 0.5))
  den <- (1 - exp(0.5 * h)) * (1 + exp(-h * (w_i - 0.5)))
  sq <- (num/den) - (g_i * x_i)
  return(sq)
}

# HIF1, NE_b2
wHIF1 <- function(NE_b2){return(NE_b2)}

#SS_b1, SS_b1
wSS_b1 <- function(SS_b1){return(SS_b1)}

#SS_b2, SS_b2
wSS_b2 <- function(SS_b2){return(SS_b2)}

#WNT7a, WNT7a
wWNT7a <- function(WNT7a){return(WNT7a)}

#NE_b1, NE_b1
wNE_b1 <- function(NE_b1){return(NE_b1)}

#NE_b2, NE_b2
wNE_b2 <- function(NE_b2){return(NE_b2)}

#Gs, PGI2 | S1P_b1 | ((NE_b1 | NE_b2) & !Superoxide)
wGs <- function(PGI2, S1P_b1, NE_b1, NE_b2, Superoxide){
  return(max(PGI2, S1P_b1, min(max(NE_b1, NE_b2), (1 - Superoxide))))
}

# AC, Gs & !Gi_b2 & !Ca_b2
wAC <- function(Gs, Gi_b2, Ca_b2){
  return(min(Gs, 1 - Gi_b2, 1 - Ca_b2))
}

# PGI2, SS_b1 & Gq & !ONOO
wPGI2 <- function(SS_b1, Gq, ONOO){
  return(min(SS_b1, Gq, 1 - ONOO))
}

# TIE2, SS_b1 & !HIF1
wTIE2 <- function(SS_b1, HIF1){
  return(min(SS_b1, 1 - HIF1))
}

# JPs, Ca_b1 & !(VEGFA_b2 & SRC_b2) & (RAC1 | WNT7a | AC)
wJPs <- function(Ca_b1, VEGFA_b2, SRC_b2, RAC1, WNT7a, AC){
  return(min(Ca_b1, 1-min(VEGFA_b2, SRC_b2), max(RAC1, WNT7a, AC)))
}

# RAC1, (S1P_b1 & (AC | NO_b1 | NO_b2) & WNT7a & SS_b1) | (!RHOA & (S1P_b1 | AC | NO_b1 | NO_b2 | WNT7a | SS_b1))
wRAC1 <- function(S1P_b1, AC, NO_b1, NO_b2, WNT7a, SS_b1, RHOA){
  a <- min(S1P_b1, max(AC, NO_b1, NO_b2), WNT7a, SS_b1)
  b <- min(1 - RHOA, max(S1P_b1, AC, NO_b1, NO_b2, WNT7a, SS_b1))
  return(max(a, b))
}

# S1P_b1, (!SRC_b2) & ETS
wS1P_b1 <- function(SRC_b2, ETS){
  return(min(1 - SRC_b2, ETS))
}

# S1P_b2, SRC_b2 & ETS
wS1P_b2 <- function(SRC_b2, ETS){
  return(min(SRC_b2, ETS))
}

# Gi_b1, S1P_b2 | S1P_b1 | NE_b2 | NE_b1
wGi_b1 <- function(S1P_b2, S1P_b1, NE_b2, NE_b1){
  return(max(S1P_b2, S1P_b1, NE_b2, NE_b1))
}

# Gi_b2, (NE_b2 & (S1P_b2 | S1P_b1)) | ( NE_b1 & Superoxide & S1P_b2 )
wGi_b2 <- function(NE_b2, S1P_b2, S1P_b1, NE_b1, Superoxide){
  return(max(min(NE_b2, max(S1P_b2, S1P_b1)), min(NE_b1, Superoxide, S1P_b2)))
}

# NO_b1, (!Ca_b2 & Ca_b1) & (!SRC_b2 | !RHOA) & (Gi_b2 | Gi_b1 | VEGFA_b2 | VEGFA_b1 | TIE2 | AC)
wNO_b1 <- function(Ca_b2, Ca_b1, SRC_b2, RHOA, Gi_b2, Gi_b1, VEGFA_b2, VEGFA_b1, TIE2, AC){
  a <- min(1 - Ca_b2, Ca_b1)
  b <- max(1 - SRC_b2, 1 - RHOA)
  c <- max(Gi_b2, Gi_b1, VEGFA_b2, VEGFA_b1, TIE2, AC)
  return(min(a, b, c))
}

# NO_b2, Ca_b2 & !SRC_b2 & !RHOA & (Gi_b2 | VEGFA_b2)
wNO_b2 <- function(Ca_b2, SRC_b2, RHOA, Gi_b2, VEGFA_b2){
  return(min(Ca_b2, 1 - SRC_b2, 1 - RHOA, max(Gi_b2, VEGFA_b2)))
}

# Ca_b1, (SS_b2 & NO_b2) | (SS_b1 & (!Gq | NO_b1 | NO_b2)) | (!SS_b2 & !SS_b1 & !NO_b2 & !Gq)
wCa_b1 <- function(SS_b2, NO_b2, SS_b1, Gq, NO_b1){ 
  a <- min(SS_b2, NO_b2)
  b <- min(SS_b1, max(1 - Gq, NO_b1, NO_b2))
  c <- min(1 - SS_b2, 1 - SS_b1, 1 - NO_b2, 1 - Gq)
  return(max(a, b, c))
}

# Ca_b2, (SS_b2 & !NO_b2) | (Gq & SS_b1 & !NO_b1 & !NO_b2)  
wCa_b2 <- function(SS_b2, NO_b2, Gq, SS_b1, NO_b1){
  return(max(min(SS_b2, 1 - NO_b2), min(Gq, SS_b1, 1 - NO_b1, 1 - NO_b2)))
}

# VEGFA_b1, (SS_b1 & !(SRC_b2 & HIF1)) | (!SS_b2 & !SS_b1 & HIF1)
wVEGFA_b1 <- function(SS_b1, SRC_b2, HIF1, SS_b2){
  a <- min(SS_b1, 1 - min (SRC_b2, HIF1))
  b <- min(1 - SS_b2, 1 - SS_b1, HIF1)
  return(max(a, b))
}

# VEGFA_b2, SS_b2 | (SS_b1 & SRC_b2 & HIF1)
wVEGFA_b2 <- function(SS_b2, SS_b1, SRC_b2, HIF1){
  return(max(SS_b2, min(SS_b1, SRC_b2, HIF1)))
}

#Gq, SS_b2 | S1P_b2 | NE_b2 | ANGII | TXA2
wGq <-function(SS_b2, S1P_b2, NE_b2, ANGII, TXA2){
  return(max(SS_b2, S1P_b2, NE_b2, ANGII, TXA2))
}

# ANGII, Gq & ETS
wANGII <- function(Gq, ETS){
  return(min(Gq, ETS))
}

# G12, ANGII | TXA2 | S1P_b2
wG12 <- function(ANGII, TXA2, S1P_b2){
  return(max(ANGII, TXA2, S1P_b2))
}

# TXA2, SS_b1 & Gq & ONOO
wTXA2 <- function(SS_b1, Gq, ONOO){
  return(min(SS_b1, Gq, ONOO))
}

# ONOO, NO_b2 & Superoxide
wONOO <- function(NO_b2, Superoxide){
  return(min(NO_b2, Superoxide))
}

# ETS, VEGFA_b2 | Gq
wETS <- function(VEGFA_b2, Gq){
  return(max(VEGFA_b2, Gq))
}

# Superoxide, !SS_b1 & (VEGFA_b2 | Gq | Ca_b2) & (HIF1 | NE_b2)
wSuperoxide <- function(SS_b1, VEGFA_b2, Gq, Ca_b2, HIF1, NE_b2){
  return(min(1 - SS_b1, max(VEGFA_b2, Gq, Ca_b2), max(HIF1, NE_b2)))
}
# SRC_b1, ((VEGFA_b2 | Gi_b2) & (TIE2 | PGI2)) | ((VEGFA_b1 | Gi_b1) & !(TIE2 | PGI2))
wSRC_b1 <- function(VEGFA_b2, Gi_b2, TIE2, PGI2, VEGFA_b1, Gi_b1){
  a <- min(max(VEGFA_b2, Gi_b2), max(TIE2, PGI2))
  b <- min(max(VEGFA_b1, Gi_b1), 1 - max(TIE2, PGI2))
  return(max(a, b))
}

# SRC_b2, (VEGFA_b2 | Gi_b2) & !TIE2 & !PGI2
wSRC_b2 <- function(VEGFA_b2, Gi_b2, TIE2, PGI2){
  return(min(max(VEGFA_b2, Gi_b2), 1 - TIE2, 1 - PGI2))
}

# RHOA, (G12 & SRC_b2) | (G12 & SRC_b1 & (!RAC1 | !AC | !WNT7a))
wRHOA <- function(G12, SRC_b2, SRC_b1, RAC1, AC, WNT7a){
  a <- min(G12, SRC_b2)
  b <- min(G12, SRC_b1, max(1 - RAC1, 1 - AC, 1- WNT7a))
  return(max(a, b))
}

# MLC, (Ca_b2 & (RHOA | !NO_b2)) | (!Ca_b2 & RHOA & !NO_b2)
wMLC <- function(Ca_b2, RHOA, NO_b2){
  a <- min(Ca_b2, max(RHOA, 1 - NO_b2))
  b <- min(1 - Ca_b2, RHOA, 1 - NO_b2)
  return(max(a, b))
}

# Define the system of ODEs
# with(as.list(c(state, parameters))
# combines the objects state and parameters into a single list, 
# then evaluates an expression within the environment of that 
# combined list. This allows access to variables from both state 
# and parameters directly by name without prefixing them.

levels_odes <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dHIF1 <- squad(HIF1, wHIF1(NE_b2), h, g)
    dSS_b1 <- squad(SS_b1, wSS_b1(SS_b1), h, g)
    dSS_b2 <- squad(SS_b2, wSS_b2(SS_b2), h, g)
    dWNT7a <- squad(WNT7a, wWNT7a(WNT7a), h, g)
    dNE_b1 <- squad(NE_b1, wNE_b1(NE_b1), h, g)
    dNE_b2 <- squad(NE_b2, wNE_b2(NE_b2), h, g)
    dGs <- squad(Gs, wGs(PGI2, S1P_b1, NE_b1, NE_b2, Superoxide), h, g)
    dAC <- squad(AC, wAC(Gs, Gi_b2, Ca_b2), h, g)
    dPGI2 <- squad(PGI2, wPGI2(SS_b1, Gq, ONOO), h, g)
    dTIE2 <- squad(TIE2, wTIE2(SS_b1, HIF1), h, g)
    dJPs <- squad(JPs, wJPs(Ca_b1, VEGFA_b2, SRC_b2, RAC1, WNT7a, AC), h, g)
    dRAC1 <- squad(RAC1, wRAC1(S1P_b1, AC, NO_b1, NO_b2, WNT7a, SS_b1, RHOA), h, g)
    dS1P_b1 <- squad(S1P_b1, wS1P_b1(SRC_b2, ETS), h, g)
    dS1P_b2 <- squad(S1P_b2, wS1P_b2(SRC_b2, ETS), h, g)
    dGi_b1 <- squad(Gi_b1, wGi_b1(S1P_b2, S1P_b1, NE_b2, NE_b1), h, g)
    dGi_b2 <- squad(Gi_b2, wGi_b2(NE_b2, S1P_b2, S1P_b1, NE_b1, Superoxide), h, g)
    dNO_b1 <- squad(NO_b1, wNO_b1(Ca_b2, Ca_b1, SRC_b2, RHOA, Gi_b2, Gi_b1, VEGFA_b2, VEGFA_b1, TIE2, AC), h, g)
    dNO_b2 <- squad(NO_b2, wNO_b2(Ca_b2, SRC_b2, RHOA, Gi_b2, VEGFA_b2), h, g)
    dCa_b1 <- squad(Ca_b1, wCa_b1(SS_b2, NO_b2, SS_b1, Gq, NO_b1), h, g)
    dCa_b2 <- squad(Ca_b2, wCa_b2(SS_b2, NO_b2, Gq, SS_b1, NO_b1), h, g)
    dVEGFA_b1 <- squad(VEGFA_b1, wVEGFA_b1(SS_b1, SRC_b2, HIF1, SS_b2), h, g)
    dVEGFA_b2 <- squad(VEGFA_b2, wVEGFA_b2(SS_b2, SS_b1, SRC_b2, HIF1), h, g)
    dGq <- squad(Gq, wGq(SS_b2, S1P_b2, NE_b2, ANGII, TXA2), h, g)
    dANGII <- squad(ANGII, wANGII(Gq, ETS), h, g)
    dG12 <- squad(G12, wG12(ANGII, TXA2, S1P_b2), h, g)
    dTXA2 <- squad(TXA2, wTXA2(SS_b1, Gq, ONOO), h, g)
    dONOO <- squad(ONOO, wONOO(NO_b2, Superoxide), h, g)
    dETS <- squad(ETS, wETS(VEGFA_b2, Gq), h, g)
    dSuperoxide <- squad(Superoxide, wSuperoxide(SS_b1, VEGFA_b2, Gq, Ca_b2, HIF1, NE_b2), h, g)
    dSRC_b1 <- squad(SRC_b1, wSRC_b1(VEGFA_b2, Gi_b2, TIE2, PGI2, VEGFA_b1, Gi_b1), h, g)
    dSRC_b2 <- squad(SRC_b2, wSRC_b2(VEGFA_b2, Gi_b2, TIE2, PGI2), h, g)
    dRHOA <- squad(RHOA, wRHOA(G12, SRC_b2, SRC_b1, RAC1, AC, WNT7a), h, g)
    dMLC <- squad(MLC, wMLC(Ca_b2, RHOA, NO_b2), h, g)
    # Return the derivatives as a list
    derivatives <- c(
      dHIF1, dSS_b1, dSS_b2,
      dWNT7a, dNE_b1, dNE_b2, 
      dGs, dAC, dPGI2, 
      dTIE2, dJPs, dRAC1, 
      dS1P_b1, dS1P_b2, dGi_b1, 
      dGi_b2, dNO_b1, dNO_b2, 
      dCa_b1, dCa_b2, dVEGFA_b1, 
      dVEGFA_b2, dGq, dANGII, 
      dG12, dTXA2, dONOO, 
      dETS, dSuperoxide, dSRC_b1, 
      dSRC_b2, dRHOA, dMLC 
      )
    return(list(derivatives))
  })
}

for(r in 1:nrow(attractors)){
  # Set the initial conditions and parameter values
  file_path <- paste('Figures\\ode\\', as.character(r), '.pdf', sep='')
  initial_state <- attractors[r, 2:ncol(attractors)]
  initial_state <- setNames(as.numeric(initial_state), names(initial_state))
  params <- c(h = 10, g = 1)
  
  # Define the time sequence to sample the solution
  times <- seq(0, 50, by = 0.001)
  
  # Solve the ODE system using the default 'lsoda' integrator
  output <- ode(y = initial_state, times = times, func = levels_odes, parms = params)

  # Save the cyclic attractor simulations as rds files
  if(r > 50){
    # Convert deSolve matrix to data frame and save it
    df_output <- as.data.frame(output)
    rds_path <- paste('Figures\\ode\\', as.character(r), '.rds', sep='')
    # Save to RDS
    saveRDS(df_output, rds_path)
  }
  
  # Open the PDF device
  pdf(file_path, width = 11, height = 8)
  
  # Generate the plots
  plot(output)
  
  # Close the device to save the file
  dev.off()
}



