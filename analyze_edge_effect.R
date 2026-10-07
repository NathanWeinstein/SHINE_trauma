## First specify the packages of interest
packages <- c("BoolNet", "ggplot2", "forcats", "cowplot")

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
if (!dir.exists("Figures\\Edges")) {
  dir.create("Figures\\Edges", recursive = TRUE)
}

make_row <- function(tss){
  tssl <- c(tss$tssv["bf"], tss$tssv_sd["bf"], 
            tss$tssv["os"], tss$tssv_sd["os"],
            tss$tssv["gd"], tss$tssv_sd["gd"],
            tss$tssv["ef"], tss$tssv_sd["ef"])
  tssl_names <- c("bf","bf_sd", "os","os_sd", "gd","gd_sd", "ef","ef_sd")
  names(tssl) <- tssl_names
  return(tssl)
}

# Analyze the wild type
initial_states <- 1000
reps <- 10

na <- analyze_net(shine, initial_states, reps)
tssv <- na$tssv
tssv_sd <- na$tssv_sd
all <- list(tssv = tssv, tssv_sd = tssv_sd)
edges <- data.frame(matrix(ncol = 8, nrow = 0))
edges <- rbind(edges, make_row(na))
rows <- c('Control')

# Analyze edges
NE_a_HIF1 <- loadNetwork("nets\\NE_a_HIF1")
a <- plot_edge_behavior(NE_a_HIF1, "NE_a_HIF1", tssv, tssv_sd)
edges <- rbind(edges, make_row(a))
rows <- c(rows, "NE_a_HIF1")

NE_a_O2 <- loadNetwork("nets\\NE_a_O2-")
b <- plot_edge_behavior(NE_a_O2, "NE_a_O2-", tssv, tssv_sd)
edges <- rbind(edges, make_row(b))
rows <- c(rows, "NE_a_O2-")

ONOO_i_PGI2 <- loadNetwork("nets\\ONOO_i_PGI2")
c <- plot_edge_behavior(ONOO_i_PGI2, "ONOO_i_PGI2", tssv, tssv_sd)
edges <- rbind(edges, make_row(c))
rows <- c(rows, "ONOO_i_PGI2")

ONOO_a_TXA2 <- loadNetwork("nets\\ONOO_a_TXA2")
d <- plot_edge_behavior(ONOO_a_TXA2, "ONOO_a_TXA2", tssv, tssv_sd)
edges <- rbind(edges, make_row(d))
rows <- c(rows, "ONOO_a_TXA2")

ETS_a_ANGII <- loadNetwork("nets\\ETS_a_ANGII")
e <- plot_edge_behavior(ETS_a_ANGII, "ETS_a_ANGII", tssv, tssv_sd)
edges <- rbind(edges, make_row(e))
rows <- c(rows, "ETS_a_ANGII")

VEGFA_a_ETS <- loadNetwork("nets\\VEGFA_a_ETS")
f <- plot_edge_behavior(VEGFA_a_ETS, "VEGFA_a_ETS", tssv, tssv_sd)
edges <- rbind(edges, make_row(f))
rows <- c(rows, "VEGFA_a_ETS")

Shear_a_TXA2 <- loadNetwork("nets\\Shear_a_TXA2")
g <- plot_edge_behavior(Shear_a_TXA2, "Shear_a_TXA2", tssv, tssv_sd)
edges <- rbind(edges, make_row(g))
rows <- c(rows, "Shear_a_TXA2")

S1P_a_G12 <- loadNetwork("nets\\S1P_a_G12")
h <- plot_edge_behavior(S1P_a_G12, "S1P_a_G12", tssv, tssv_sd)
edges <- rbind(edges, make_row(h))
rows <- c(rows, "S1P_a_G12")

Shear_i_O2 <- loadNetwork("nets\\Shear_i_O2-")
i <- plot_edge_behavior(Shear_i_O2, "Shear_i_O2-", tssv, tssv_sd)
edges <- rbind(edges, make_row(i))
rows <- c(rows, "Shear_i_O2-")

ANGII_a_G12 <- loadNetwork("nets\\ANGII_a_G12")
j <- plot_edge_behavior(ANGII_a_G12, "ANGII_a_G12", tssv, tssv_sd)
edges <- rbind(edges, make_row(j))
rows <- c(rows, "ANGII_a_G12")

TXA2_a_G12 <- loadNetwork("nets\\TXA2_a_G12")
k <- plot_edge_behavior(TXA2_a_G12, "TXA2_a_G12", tssv, tssv_sd)
edges <- rbind(edges, make_row(k))
rows <- c(rows, "TXA2_a_G12")

O2_a_ONOO <- loadNetwork("nets\\O2-_a_ONOO")
l <- plot_edge_behavior(O2_a_ONOO, "O2-_a_ONOO", tssv, tssv_sd)
edges <- rbind(edges, make_row(l))
rows <- c(rows, "O2-_a_ONOO")

HIF1_a_O2 <- loadNetwork("nets\\HIF1_a_O2-")
m <- plot_edge_behavior(HIF1_a_O2, "HIF1_a_O2-", tssv, tssv_sd)
edges <- rbind(edges, make_row(m))
rows <- c(rows, "HIF1_a_O2-")

TXA2_a_Gq <- loadNetwork("nets\\TXA2_a_Gq")
n <- plot_edge_behavior(TXA2_a_Gq, "TXA2_a_Gq", tssv, tssv_sd)
edges <- rbind(edges, make_row(n))
rows <- c(rows, "TXA2_a_Gq")

Gq_a_TXA2 <- loadNetwork("nets\\Gq_a_TXA2")
o <- plot_edge_behavior(Gq_a_TXA2, "Gq_a_TXA2", tssv, tssv_sd)
edges <- rbind(edges, make_row(o))
rows <- c(rows, "Gq_a_TXA2")

Gq_a_ETS <- loadNetwork("nets\\Gq_a_ETS")
p <- plot_edge_behavior(Gq_a_ETS, "Gq_a_ETS", tssv, tssv_sd)
edges <- rbind(edges, make_row(p))
rows <- c(rows, "Gq_a_ETS")

Gq_a_O2 <- loadNetwork("nets\\Gq_a_O2-")
q <- plot_edge_behavior(Gq_a_O2, "Gq_a_O2-", tssv, tssv_sd)
edges <- rbind(edges, make_row(q))
rows <- c(rows, "Gq_a_O2-")

Shear_a_Gq <- loadNetwork("nets\\Shear_a_Gq")
r <- plot_edge_behavior(Shear_a_Gq, "Shear_a_Gq", tssv, tssv_sd)
edges <- rbind(edges, make_row(r))
rows <- c(rows, "Shear_a_Gq")

VEGFA_a_SRC <- loadNetwork("nets\\VEGFA_a_SRC")
s <- plot_edge_behavior(VEGFA_a_SRC, "VEGFA_a_SRC", tssv, tssv_sd)
edges <- rbind(edges, make_row(s))
rows <- c(rows, "VEGFA_a_SRC")

G12_a_RHOA <- loadNetwork("nets\\G12_a_RHOA")
t <- plot_edge_behavior(G12_a_RHOA, "G12_a_RHOA", tssv, tssv_sd)
edges <- rbind(edges, make_row(t))
rows <- c(rows, "G12_a_RHOA")

SRC_a_RHOA <- loadNetwork("nets\\SRC_a_RHOA")
u <- plot_edge_behavior(SRC_a_RHOA, "SRC_a_RHOA", tssv, tssv_sd)
edges <- rbind(edges, make_row(u))
rows <- c(rows, "SRC_a_RHOA")

NO_a_ONOO <- loadNetwork("nets\\NO_a_ONOO")
v <- plot_edge_behavior(NO_a_ONOO, "NO_a_ONOO", tssv, tssv_sd)
edges <- rbind(edges, make_row(v))
rows <- c(rows, "NO_a_ONOO")

WNT7_i_RHOA <- loadNetwork("nets\\WNT7_i_RHOA")
w <- plot_edge_behavior(WNT7_i_RHOA, "WNT7_i_RHOA", tssv, tssv_sd)
edges <- rbind(edges, make_row(w))
rows <- c(rows, "WNT7_i_RHOA")

Ca2_a_O2 <- loadNetwork("nets\\Ca2_a_O2")
x <- plot_edge_behavior(Ca2_a_O2, "Ca2_a_O2-", tssv, tssv_sd)
edges <- rbind(edges, make_row(x))
rows <- c(rows, "Ca2_a_O2-")

VEGFA_a_O2 <- loadNetwork("nets\\VEGFA_a_O2")
y <- plot_edge_behavior(VEGFA_a_O2, "VEGFA_a_O2-", tssv, tssv_sd)
edges <- rbind(edges, make_row(y))
rows <- c(rows, "VEGFA_a_O2-")

Gq_a_ANGII <- loadNetwork("nets\\Gq_a_ANGII")
z <- plot_edge_behavior(Gq_a_ANGII, "Gq_a_ANGII", tssv, tssv_sd)
edges <- rbind(edges, make_row(z))
rows <- c(rows, "Gq_a_ANGII")

ANGII_a_Gq <- loadNetwork("nets\\ANGII_a_Gq")
aa <- plot_edge_behavior(ANGII_a_Gq, "ANGII_a_Gq", tssv, tssv_sd)
edges <- rbind(edges, make_row(aa))
rows <- c(rows, "ANGII_a_Gq")

SRC_i_JPs <- loadNetwork("nets\\SRC_i_JPs")
ab <- plot_edge_behavior(SRC_i_JPs, "SRC_i_JPs", tssv, tssv_sd)
edges <- rbind(edges, make_row(ab))
rows <- c(rows, "SRC_i_JPs")

AC_i_RHOA <- loadNetwork("nets\\AC_i_RHOA")
ac <- plot_edge_behavior(AC_i_RHOA, "AC_i_RHOA", tssv, tssv_sd)
edges <- rbind(edges, make_row(ac))
rows <- c(rows, "AC_i_RHOA")

RHOA_a_MLC <- loadNetwork("nets\\RHOA_a_MLC")
ad <- plot_edge_behavior(RHOA_a_MLC, "RHOA_a_MLC", tssv, tssv_sd)
edges <- rbind(edges, make_row(ad))
rows <- c(rows, "RHOA_a_MLC")

HIF1_i_TIE2 <- loadNetwork("nets\\HIF1_i_TIE2")
ae <- plot_edge_behavior(HIF1_i_TIE2, "HIF1_i_TIE2", tssv, tssv_sd)
edges <- rbind(edges, make_row(ae))
rows <- c(rows, "HIF1_i_TIE2")

S1P_a_Gq <- loadNetwork("nets\\S1P_a_Gq")
af <- plot_edge_behavior(S1P_a_Gq, "S1P_a_Gq", tssv, tssv_sd)
edges <- rbind(edges, make_row(af))
rows <- c(rows, "S1P_a_Gq")

Gq_a_PGI2 <- loadNetwork("nets\\Gq_a_PGI2")
ag <- plot_edge_behavior(Gq_a_PGI2, "Gq_a_PGI2", tssv, tssv_sd)
edges <- rbind(edges, make_row(ag))
rows <- c(rows, "Gq_a_PGI2")

SRC_a_S1P <- loadNetwork("nets\\SRC_a_S1P")
ah <- plot_edge_behavior(SRC_a_S1P, "SRC_a_S1P", tssv, tssv_sd)
edges <- rbind(edges, make_row(ah))
rows <- c(rows, "SRC_a_S1P")

NE_a_Gq <- loadNetwork("nets\\NE_a_Gq")
ai <- plot_edge_behavior(NE_a_Gq, "NE_a_Gq", tssv, tssv_sd)
edges <- rbind(edges, make_row(ai))
rows <- c(rows, "NE_a_Gq")

Gq_a_Ca2 <- loadNetwork("nets\\Gq_a_Ca2")
aj <- plot_edge_behavior(Gq_a_Ca2, "Gq_a_Ca2", tssv, tssv_sd)
edges <- rbind(edges, make_row(aj))
rows <- c(rows, "Gq_a_Ca2")

Gi_a_SRC <- loadNetwork("nets\\Gi_a_SRC")
ak <- plot_edge_behavior(Gi_a_SRC, "Gi_a_SRC", tssv, tssv_sd)
edges <- rbind(edges, make_row(ak))
rows <- c(rows, "Gi_a_SRC")

PGI2_i_SRC <- loadNetwork("nets\\PGI2_i_SRC")
al <- plot_edge_behavior(PGI2_i_SRC, "PGI2_i_SRC", tssv, tssv_sd)
edges <- rbind(edges, make_row(al))
rows <- c(rows, "PGI2_i_SRC")

TIE2_i_SRC <- loadNetwork("nets\\TIE2_i_SRC")
am <- plot_edge_behavior(TIE2_i_SRC, "TIE2_i_SRC", tssv, tssv_sd)
edges <- rbind(edges, make_row(am))
rows <- c(rows, "TIE2_i_SRC")

SRC_i_NO <- loadNetwork("nets\\SRC_i_NO")
an <- plot_edge_behavior(SRC_i_NO, "SRC_i_NO", tssv, tssv_sd)
edges <- rbind(edges, make_row(an))
rows <- c(rows, "SRC_i_NO")

RHOA_i_RAC1 <- loadNetwork("nets\\RHOA_i_RAC1")
ao <- plot_edge_behavior(RHOA_i_RAC1, "RHOA_i_RAC1", tssv, tssv_sd)
edges <- rbind(edges, make_row(ao))
rows <- c(rows, "RHOA_i_RAC1")

RHOA_i_NO <- loadNetwork("nets\\RHOA_i_NO")
ap <- plot_edge_behavior(RHOA_i_NO, "RHOA_i_NO", tssv, tssv_sd)
edges <- rbind(edges, make_row(ap))
rows <- c(rows, "RHOA_i_NO")

RAC1_i_RHOA <- loadNetwork("nets\\RAC1_i_RHOA")
aq <- plot_edge_behavior(RAC1_i_RHOA, "RAC1_i_RHOA", tssv, tssv_sd)
edges <- rbind(edges, make_row(aq))
rows <- c(rows, "RAC1_i_RHOA")

NO_i_MLC <- loadNetwork("nets\\NO_i_MLC")
ar <- plot_edge_behavior(NO_i_MLC, "NO_i_MLC", tssv, tssv_sd)
edges <- rbind(edges, make_row(ar))
rows <- c(rows, "NO_i_MLC")

HIF1_a_VEGFA <- loadNetwork("nets\\HIF1_a_VEGFA")
as <- plot_edge_behavior(HIF1_a_VEGFA, "HIF1_a_VEGFA", tssv, tssv_sd)
edges <- rbind(edges, make_row(as))
rows <- c(rows, "HIF1_a_VEGFA")

S1P_a_RAC1 <- loadNetwork("nets\\S1P_a_RAC1")
at <- plot_edge_behavior(S1P_a_RAC1, "S1P_a_RAC1", tssv, tssv_sd)
edges <- rbind(edges, make_row(at))
rows <- c(rows, "S1P_a_RAC1")

S1P_a_Gs <- loadNetwork("nets\\S1P_a_Gs")
au <- plot_edge_behavior(S1P_a_Gs, "S1P_a_Gs", tssv, tssv_sd)
edges <- rbind(edges, make_row(au))
rows <- c(rows, "S1P_a_Gs")

S1P_a_Gi <- loadNetwork("nets\\S1P_a_Gi")
av <- plot_edge_behavior(S1P_a_Gi, "S1P_a_Gi", tssv, tssv_sd)
edges <- rbind(edges, make_row(av))
rows <- c(rows, "S1P_a_Gi")

NE_a_Gs <- loadNetwork("nets\\NE_a_Gs")
aw <- plot_edge_behavior(NE_a_Gs, "NE_a_Gs", tssv, tssv_sd)
edges <- rbind(edges, make_row(aw))
rows <- c(rows, "NE_a_Gs")

NE_a_Gi <- loadNetwork("nets\\NE_a_Gi")
ax <- plot_edge_behavior(NE_a_Gi, "NE_a_Gi", tssv, tssv_sd)
edges <- rbind(edges, make_row(ax))
rows <- c(rows, "NE_a_Gi")

Gi_a_NO <- loadNetwork("nets\\Gi_a_NO")
ay <- plot_edge_behavior(Gi_a_NO, "Gi_a_NO", tssv, tssv_sd)
edges <- rbind(edges, make_row(ay))
rows <- c(rows, "Gi_a_NO")

Ca_a_NO <- loadNetwork("nets\\Ca_a_NO")
az <- plot_edge_behavior(Ca_a_NO, "Ca_a_NO", tssv, tssv_sd)
edges <- rbind(edges, make_row(az))
rows <- c(rows, "Ca_a_NO")

NO_i_Ca <- loadNetwork("nets\\NO_i_Ca")
ba <- plot_edge_behavior(NO_i_Ca, "NO_i_Ca", tssv, tssv_sd)
edges <- rbind(edges, make_row(ba))
rows <- c(rows, "NO_i_Ca")

TIE2_a_NO <- loadNetwork("nets\\TIE2_a_NO")
bb <- plot_edge_behavior(TIE2_a_NO, "TIE2_a_NO", tssv, tssv_sd)
edges <- rbind(edges, make_row(bb))
rows <- c(rows, "TIE2_a_NO")

Gi_i_AC <- loadNetwork("nets\\Gi_i_AC")
bc <- plot_edge_behavior(Gi_i_AC, "Gi_i_AC", tssv, tssv_sd)
edges <- rbind(edges, make_row(bc))
rows <- c(rows, "Gi_i_AC")

AC_a_NO <- loadNetwork("nets\\AC_a_NO")
bd <- plot_edge_behavior(AC_a_NO, "AC_a_NO", tssv, tssv_sd)
edges <- rbind(edges, make_row(bd))
rows <- c(rows, "AC_a_NO")

NO_a_RAC1 <- loadNetwork("nets\\NO_a_RAC1")
be <- plot_edge_behavior(NO_a_RAC1, "NO_a_RAC1", tssv, tssv_sd)
edges <- rbind(edges, make_row(be))
rows <- c(rows, "NO_a_RAC1")

VEGFA_a_NO <- loadNetwork("nets\\VEGFA_a_NO")
bf <- plot_edge_behavior(VEGFA_a_NO, "VEGFA_a_NO", tssv, tssv_sd)
edges <- rbind(edges, make_row(bf))
rows <- c(rows, "VEGFA_a_NO")

Ca_i_AC <- loadNetwork("nets\\Ca_i_AC")
bg <- plot_edge_behavior(Ca_i_AC, "Ca_i_AC", tssv, tssv_sd)
edges <- rbind(edges, make_row(bg))
rows <- c(rows, "Ca_i_AC")

Ca_a_AC <- loadNetwork("nets\\Ca_a_AC_e")
bg_e <- plot_edge_behavior(Ca_a_AC, "Ca_a_AC8", tssv, tssv_sd)
edges <- rbind(edges, make_row(bg_e))
rows <- c(rows, "Ca_a_AC8")

Ca_a_JPs <- loadNetwork("nets\\Ca_a_JPs")
bh <- plot_edge_behavior(Ca_a_JPs, "Ca_a_JPs", tssv, tssv_sd)
edges <- rbind(edges, make_row(bh))
rows <- c(rows, "Ca_a_JPs")

Shear_a_Ca <- loadNetwork("nets\\Shear_a_Ca")
bi <- plot_edge_behavior(Shear_a_Ca, "Shear_a_Ca", tssv, tssv_sd)
edges <- rbind(edges, make_row(bi))
rows <- c(rows, "Shear_a_Ca")

VEGFA_i_JPs <- loadNetwork("nets\\VEGFA_i_JPs")
bj <- plot_edge_behavior(VEGFA_i_JPs, "VEGFA_i_JPs", tssv, tssv_sd)
edges <- rbind(edges, make_row(bj))
rows <- c(rows, "VEGFA_i_JPs")

Shear_a_VEGFA <- loadNetwork("nets\\Shear_a_VEGFA")
bk <- plot_edge_behavior(Shear_a_VEGFA, "Shear_a_VEGFA", tssv, tssv_sd)
edges <- rbind(edges, make_row(bk))
rows <- c(rows, "Shear_a_VEGFA")

Shear_a_TIE2 <- loadNetwork("nets\\Shear_a_TIE2")
bl <- plot_edge_behavior(Shear_a_TIE2, "Shear_a_TIE2", tssv, tssv_sd)
edges <- rbind(edges, make_row(bl))
rows <- c(rows, "Shear_a_TIE2")

Shear_a_PGI2 <- loadNetwork("nets\\Shear_a_PGI2")
bm <- plot_edge_behavior(Shear_a_PGI2, "Shear_a_PGI2", tssv, tssv_sd)
edges <- rbind(edges, make_row(bm))
rows <- c(rows, "Shear_a_PGI2")

PGI2_a_Gs <- loadNetwork("nets\\PGI2_a_Gs")
bn <- plot_edge_behavior(PGI2_a_Gs, "PGI2_a_Gs", tssv, tssv_sd)
edges <- rbind(edges, make_row(bn))
rows <- c(rows, "PGI2_a_Gs")

Gs_a_AC <- loadNetwork("nets\\Gs_a_AC")
bo <- plot_edge_behavior(Gs_a_AC, "Gs_a_AC", tssv, tssv_sd)
edges <- rbind(edges, make_row(bo))
rows <- c(rows, "Gs_a_AC")

AC_a_RAC1 <- loadNetwork("nets\\AC_a_RAC1")
bp <- plot_edge_behavior(AC_a_RAC1, "AC_a_RAC1", tssv, tssv_sd)
edges <- rbind(edges, make_row(bp))
rows <- c(rows, "AC_a_RAC1")

AC_a_JPs <- loadNetwork("nets\\AC_a_JPs")
bq <- plot_edge_behavior(AC_a_JPs, "AC_a_JPs", tssv, tssv_sd)
edges <- rbind(edges, make_row(bq))
rows <- c(rows, "AC_a_JPs")

Shear_a_RAC1 <- loadNetwork("nets\\Shear_a_RAC1")
br <- plot_edge_behavior(Shear_a_RAC1, "Shear_a_RAC1", tssv, tssv_sd)
edges <- rbind(edges, make_row(br))
rows <- c(rows, "Shear_a_RAC1")

RAC1_a_JPs <- loadNetwork("nets\\RAC1_a_JPs")
bs <- plot_edge_behavior(RAC1_a_JPs, "RAC1_a_JPs", tssv, tssv_sd)
edges <- rbind(edges, make_row(bs))
rows <- c(rows, "RAC1_a_JPs")

WNT7a_a_RAC1 <- loadNetwork("nets\\WNT7a_a_RAC1")
bt <- plot_edge_behavior(WNT7a_a_RAC1, "WNT7a_a_RAC1", tssv, tssv_sd)
edges <- rbind(edges, make_row(bt))
rows <- c(rows, "WNT7a_a_RAC1")

WNT7a_a_JPs <- loadNetwork("nets\\WNT7a_a_JPs")
bu <- plot_edge_behavior(WNT7a_a_JPs, "WNT7a_a_JPs", tssv, tssv_sd)
edges <- rbind(edges, make_row(bu))
rows <- c(rows, "WNT7a_a_JPs")

Superoxide_i_Gs <- loadNetwork("nets\\Superoxide_i_Gs")
bv <- plot_edge_behavior(Superoxide_i_Gs, "Superoxide_i_Gs", tssv, tssv_sd)
edges <- rbind(edges, make_row(bv))
rows <- c(rows, "Superoxide_i_Gs")

ETS_a_S1P <- loadNetwork("nets\\ETS_a_S1P")
bw <- plot_edge_behavior(ETS_a_S1P, "ETS_a_S1P", tssv, tssv_sd)
edges <- rbind(edges, make_row(bw))
rows <- c(rows, "ETS_a_S1P")

Superoxide_a_Gi <- loadNetwork("nets\\Superoxide_a_Gi")
bx <- plot_edge_behavior(Superoxide_a_Gi, "Superoxide_a_Gi", tssv, tssv_sd)
edges <- rbind(edges, make_row(bx))
rows <- c(rows, "Superoxide_a_Gi")

SRC_a_VEGFA <- loadNetwork("nets\\SRC_a_VEGFA")
by <- plot_edge_behavior(SRC_a_VEGFA, "SRC_a_VEGFA", tssv, tssv_sd)
edges <- rbind(edges, make_row(by))
rows <- c(rows, "SRC_a_VEGFA")

Ca_a_MLC <- loadNetwork("nets\\Ca_a_MLC")
bz <- plot_edge_behavior(Ca_a_MLC, "Ca_a_MLC", tssv, tssv_sd)
edges <- rbind(edges, make_row(bz))
rows <- c(rows, "Ca_a_MLC")

colnames(edges) <- c("bf","bf_sd", "os","os_sd", "gd","gd_sd", "ef","ef_sd")
rownames(edges) <- rows
write.csv(edges, "edges.csv")

# Barrier function
edges_bf <- edges[order(edges$bf), ]
a <- ggplot(edges_bf) +
  geom_bar(aes(x=fct_inorder(rownames(edges_bf)), y=bf, fill=bf), stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rownames(edges_bf), ymin=bf-bf_sd, ymax=bf+bf_sd), width=0.4, colour="orange", alpha=0.9, linewidth=1.3) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Edge removed", y = "Barrier function")

# Oxidative stress
edges_os <- edges[order(edges$os), ]
b <- ggplot(edges_os) +
  geom_bar(aes(x=fct_inorder(rownames(edges_os)), y=os, fill=os), stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rownames(edges_os), ymin=os-os_sd, ymax=os+os_sd), width=0.4, colour="orange", alpha=0.9, linewidth=1.3) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Edge removed", y = "Oxidative stress")

# Glycocalyx damage
edges_gd <- edges[order(edges$gd), ]
c <- ggplot(edges_gd) +
  geom_bar(aes(x=fct_inorder(rownames(edges_gd)), y=gd, fill=gd), stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rownames(edges_gd), ymin=gd-gd_sd, ymax=gd+gd_sd), width=0.4, colour="orange", alpha=0.9, linewidth=1.3) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Edge removed", y = "Oxidative stress")

# Endothelial function
edges_ef <- edges[order(edges$ef), ]
d <- ggplot(edges_ef) +
  geom_bar(aes(x=fct_inorder(rownames(edges_ef)), y=ef, fill=ef), stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rownames(edges_ef), ymin=ef-ef_sd, ymax=ef+ef_sd), width=0.4, colour="orange", alpha=0.9, linewidth=1.3) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Edge removed", y = "Endothelial function")

jpeg(filename = "Figures\\main\\edge_results.jpeg", 
     width = 1800, height = 1800, quality = 95)
plot_grid(a, b, c, d,  
          labels = c("A", "B", "C", "D"),
          ncol = 1, nrow = 4)
dev.off()