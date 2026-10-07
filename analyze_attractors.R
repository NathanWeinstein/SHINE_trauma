## First specify the packages of interest
packages <- c("BoolNet", "gt", "graphics", "grDevices", "eulerr")

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

attr_info <- getAttractors(
  shine,
  type = "synchronous",
  method = "sat.exhaustive"
)


atr1 <- getAttractorSequence(attr_info, 1)
nc <- dim(atr1)[2] + 1

attractors <- data.frame(matrix(ncol = nc, nrow = 0))
for(i in 1:length(attr_info$attractors)){
  atri <- getAttractorSequence(attr_info, i)
  atri <- cbind(attractor = rep(i, times = nrow(atri)), atri)
  attractors <- rbind(attractors, atri)
}

if (!dir.exists("Figures")) {
  dir.create("Figures", recursive = TRUE)
}
if (!dir.exists("Figures\\main")) {
  dir.create("Figures\\main", recursive = TRUE)
}
filepath <- "Figures\\main\\attractors.html"
save_attractors <- function(attractors, filepath)
  {
    gttable <- gt(attractors)
    gttable <- tab_options(
      data = gttable,
      table.font.size = px(14),
      data_row.padding = px(0.1),
      column_labels.padding = px(0)
    )
    gttable <- data_color(
      data = gttable,
      columns = 2:last_col(), 
      method = "bin", 
      bins = 2,
      palette = c("#447dae", "#f5e9a6"),
      domain = c(0, 1))
    gtsave(gttable, filepath)
    return(gttable)
  }

gttable <- save_attractors(attractors, filepath) 

make_attractor_heatmap <- function(attractors)
{
  x  <- as.matrix(attractors[,2:ncol(attractors)])
  pdf("Figures\\main\\heatmap.pdf", width = 8, height = 11)
  hv <- heatmap(x,
                col = terrain.colors(256), 
                scale = "column",
                margins = c(5,5),
                xlab = "Molecule",
                cexRow = 0.4,
                ylab = "State",
                cexCol = 0.4)
  dev.off()
  return(hv)
}

hv <- make_attractor_heatmap(attractors)

behaviors <- analyze_attractors(shine, attr_info)
bf <- behaviors[behaviors$bf == 1, "i"]
os <- behaviors[behaviors$os == 1, "i"]
gd <- behaviors[behaviors$gd == 1, "i"]
ef <- behaviors[behaviors$ef == 1, "i"]
sbf <- sum_attractors(bf, shine, attr_info)
sef <- sum_attractors(ef, shine, attr_info)
sgd <- sum_attractors(gd, shine, attr_info)

jps <- unique(attractors[attractors$JPs == 1, "attractor"])
pef_select <- attractors$NO_b1 == 1 | attractors$NO_b2 == 1 | attractors$PGI2 == 1
pef <- unique(attractors[pef_select, "attractor"])


# use list as input
A <-list(bf = bf,
         ef = ef,
         os = os,
         gd = gd,
         jps = jps,
         pef = pef)

# create venn diagram and display all the sets
set.seed(1)
print(plot(euler(A, shape = "ellipse"), quantities = TRUE)) 

save(shine, 
     attr_info, 
     attractors, 
     behaviors,
     ef,
     os,
     gd,
     bf,
     pef,
     jps,
     file = "trap_spaces.RData")


