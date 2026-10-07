## First specify the packages of interest
packages <- c("BoolNet")

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
shine <- loadNetwork("nets\\levels.bnet")
loaded <- load("all_attractors.RData")
specs <- c("JPs" = 1)

# Data frame to store the analysis
behaviors <- data.frame(matrix(ncol = 18, nrow = 0))

fields <- c(
  "i",
  "order",
  "NE", 
  "SS",
  "HIF1", 
  "WNT7a",
  "Ca_b1",
  "Ca_b2",
  "NO_b1",
  "NO_b2",
  "PGI2",
  "TXA2",
  "jps",
  "mlc",
  "bf", 
  "os", 
  "gd", 
  "ef"
)



for(i in 1:length(attr_info$attractors)){
  atri <- getAttractorSequence(attr_info, i)
  order <- nrow(atri)
  if(atri[, "NE_b2"][1] == 1){
    ne <- 2
  }else if(atri[, "NE_b1"][1] == 1){
    ne <- 1
  }else{
    ne <- 0
  }
  if(atri[, "SS_b2"][1] == 1){
    ss <- 2
  }else if(atri[, "SS_b1"][1] == 1){
    ss <- 1
  }else{
    ss <- 0
  }
  hif1 <- atri[, "HIF1"][1]
  wnt7a <- atri[, "WNT7a"][1]
  ca_b1 <- verify_attractor(atri, c(Ca_b1 = 1))
  ca_b2 <- verify_attractor(atri, c(Ca_b2 = 1))
  jps <- verify_attractor(atri, c(JPs = 1))
  mlc <- verify_attractor(atri, c(MLC=0))
  bf <- verify_attractor(atri, c(JPs = 1, MLC=0))
  os <- verify_attractor(atri, c(Superoxide=1))
  e <- verify_attractor(atri, c(ETS=1))
  gd <- os && e
  no_b1 <- verify_attractor(atri, c(NO_b1 = 1))
  no_b2 <- verify_attractor(atri, c(NO_b2 = 1))
  pgi2 <- verify_attractor(atri, c(PGI2=1))
  txa2 <- verify_attractor(atri, c(TXA2=1))
  ef <- (!txa2 && (pgi2 || no_b1 || no_b2)) || (pgi2 && (no_b1 || no_b2))
  obs <- c(
    "i" = i,
    "order" = order,
    "Cat" = ne, 
    "SS" = ss,
    "HIF1" = hif1, 
    "WNT7a" = wnt7a,
    "Ca_b1" = ca_b1,
    "Ca_b2" = ca_b2,
    "NO_b1" = no_b1,
    "NO_b2" = no_b2,
    "PGI2" = pgi2,
    "TXA2" = txa2,
    "jps" = jps,
    "mlc" = mlc,
    "bf" = bf, 
    "os" = os, 
    "gd" = gd, 
    "ef" = ef
  )
  behaviors <- rbind(behaviors, obs)
}

colnames(behaviors) <- fields
write.csv(behaviors, "levels_analysis.csv")