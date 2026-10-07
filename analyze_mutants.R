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
initial_states <- 1000
repetitions <- 10
tss <- analyze_net(shine, initial_states, repetitions)
tssv <- c(tss$tssv["bf"], tss$tssv_sd["bf"], 
          tss$tssv["os"], tss$tssv_sd["os"],
          tss$tssv["gd"], tss$tssv_sd["gd"],
          tss$tssv["ef"], tss$tssv_sd["ef"])
tssv_names <- c("bf","bf_sd", "os","os_sd", "gd","gd_sd", "ef","ef_sd")
names(tssv) <- tssv_names

# Explore loss of function mutations
lf_tss <- data.frame(matrix(ncol = 8, nrow = 0))
lf_tss <- rbind(lf_tss, tssv)
for(g in shine$genes){
  m_shine <- fixGenes(shine, g, 0)
  m_tss <- analyze_net(m_shine, initial_states, repetitions)
  m_tssv <- c(m_tss$tssv["bf"], m_tss$tssv_sd["bf"], 
              m_tss$tssv["os"], m_tss$tssv_sd["os"],
              m_tss$tssv["gd"], m_tss$tssv_sd["gd"],
              m_tss$tssv["ef"], m_tss$tssv_sd["ef"])
  lf_tss <- rbind(lf_tss, m_tssv)
}
colnames(lf_tss) <- names(tssv)
rownames(lf_tss) <- c("Control", shine$genes)
write.csv(lf_tss, "lf_tss.csv")

# Barrier function
lf_bf <- lf_tss[order(lf_tss$bf), ]
a <- ggplot(lf_bf) +
  geom_bar(aes(x=fct_inorder(rownames(lf_bf)), y=bf, fill=bf), stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rownames(lf_bf), ymin=bf-bf_sd, ymax=bf+bf_sd), width=0.4, colour="orange", alpha=0.9, linewidth=1.3) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Loss of function", y = "Barrier function")

# Oxidative stress
lf_os <- lf_tss[order(lf_tss$os), ]
c <- ggplot(lf_os) +
  geom_bar(aes(x=fct_inorder(rownames(lf_os)), y=os, fill=os), stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rownames(lf_os), ymin=os-os_sd, ymax=os+os_sd), width=0.4, colour="orange", alpha=0.9, linewidth=1.3) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Loss of function", y = "Oxidative stress")

# Glycocalyx damage
lf_gd <- lf_tss[order(lf_tss$gd), ]
e <- ggplot(lf_gd) +
  geom_bar(aes(x=fct_inorder(rownames(lf_gd)), y=gd, fill=gd), stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rownames(lf_gd), ymin=gd-gd_sd, ymax=gd+gd_sd), width=0.4, colour="orange", alpha=0.9, linewidth=1.3)+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Loss of function", y = "Glycocalyx damage")

# Endothelial function
lf_ef <- lf_tss[order(lf_tss$ef), ]
g <- ggplot(lf_ef) +
  geom_bar(aes(x=fct_inorder(rownames(lf_ef)), y=ef, fill=ef), stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rownames(lf_ef), ymin=ef-ef_sd, ymax=ef+ef_sd), width=0.4, colour="orange", alpha=0.9, linewidth=1.3)+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Loss of function", y = "Vasodilator production")


# Explore gain of function mutations
gf_tss <- data.frame(matrix(ncol = 8, nrow = 0))
gf_tss <- rbind(gf_tss, tssv)
for(g in shine$genes){
  m_shine <- fixGenes(shine, g, 1)
  m_tss <- analyze_net(m_shine, initial_states, repetitions)
  m_tssv <- c(m_tss$tssv["bf"], m_tss$tssv_sd["bf"], 
              m_tss$tssv["os"], m_tss$tssv_sd["os"],
              m_tss$tssv["gd"], m_tss$tssv_sd["gd"],
              m_tss$tssv["ef"], m_tss$tssv_sd["ef"])
  gf_tss <- rbind(gf_tss, m_tssv)
}
colnames(gf_tss) <- names(tssv)
rownames(gf_tss) <- c("Control", shine$genes)
write.csv(gf_tss, "gf_tss.csv")

# Barrier function
gf_bf <- gf_tss[order(gf_tss$bf), ]
b <- ggplot(gf_bf) +
  geom_bar(aes(x=fct_inorder(rownames(gf_bf)), y=bf, fill=bf), stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rownames(gf_bf), ymin=bf-bf_sd, ymax=bf+bf_sd), width=0.4, colour="orange", alpha=0.9, linewidth=1.3)+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Constitutive activity", y = "Barrier function")

# Oxidative stress
gf_os <- gf_tss[order(gf_tss$os), ]
d <- ggplot(gf_os) +
  geom_bar(aes(x=fct_inorder(rownames(gf_os)), y=os, fill=os), stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rownames(gf_os), ymin=os-os_sd, ymax=os+os_sd), width=0.4, colour="orange", alpha=0.9, linewidth=1.3)+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Constitutive activity", y = "Oxidative stress")

# Glycocalyx damage
gf_gd <- gf_tss[order(gf_tss$gd), ]
f <- ggplot(gf_gd) +
  geom_bar(aes(x=fct_inorder(rownames(gf_gd)), y=gd, fill=gd), stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rownames(gf_gd), ymin=gd-gd_sd, ymax=gd+gd_sd), width=0.4, colour="orange", alpha=0.9, linewidth=1.3)+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Constitutive activity", y = "Glycocalyx damage")


# Endothelial function
gf_ef <- gf_tss[order(gf_tss$ef), ]
h <- ggplot(gf_ef) +
  geom_bar(aes(x=fct_inorder(rownames(gf_ef)), y=ef, fill=ef), stat="identity", alpha=0.7) +
  geom_errorbar(aes(x=rownames(gf_ef), ymin=ef-ef_sd, ymax=ef+ef_sd), width=0.4, colour="orange", alpha=0.9, linewidth=1.3)+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(x = "Constitutive activity", y = "Vasodilator production")

jpeg(filename = "Figures\\main\\gf_lf_results.jpeg", 
     width = 1200, height = 1800, quality = 95)
plot_grid(a, b, c, d, e, f, g, h, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 1, nrow = 8)
dev.off()


