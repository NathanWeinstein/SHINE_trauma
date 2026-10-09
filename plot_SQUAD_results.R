## First specify the packages of interest
packages <- c("ggplot2", "patchwork")

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

df51 <- readRDS("Figures\\ode\\51.rds")
df52 <- readRDS("Figures\\ode\\52.rds")

plot51_AC <- ggplot() +
  geom_line(data = df51, aes(x = 1:nrow(df51), y = AC, color = "55a")) +
  geom_line(data = df52, aes(x = 1:nrow(df52), y = AC, color = "55b")) +
  scale_color_manual(values = c("55a" = "blue", "55b" = "green")) +
  theme_minimal() +
  labs(x = "time", y = "AC activity", color = "Attractor state")

plot51_JPs <- ggplot() +
  geom_line(data = df51, aes(x = 1:nrow(df51), y = JPs, color = "55a")) +
  geom_line(data = df52, aes(x = 1:nrow(df52), y = JPs, color = "55b")) +
  scale_color_manual(values = c("55a" = "blue", "55b" = "green")) +
  theme_minimal() +
  labs(x = "time", y = "JPs function", color = "Attractor state")

plot51_NO <- ggplot() +
  geom_line(data = df51, aes(x = 1:nrow(df51), y = NO_b1, color = "55a")) +
  geom_line(data = df52, aes(x = 1:nrow(df52), y = NO_b1, color = "55b")) +
  scale_color_manual(values = c("55a" = "blue", "55b" = "green")) +
  theme_minimal() +
  labs(x = "time", y = "NO production", color = "Attractor state")

plot51_Ca1 <- ggplot() +
  geom_line(data = df51, aes(x = 1:nrow(df51), y = Ca_b1, color = "55a")) +
  geom_line(data = df52, aes(x = 1:nrow(df52), y = Ca_b1, color = "55b")) +
  scale_color_manual(values = c("55a" = "blue", "55b" = "green")) +
  theme_minimal() +
  labs(x = "time", y = "Normal calcium", color = "Attractor state")

plot51_Ca2 <- ggplot() +
  geom_line(data = df51, aes(x = 1:nrow(df51), y = Ca_b2, color = "55a")) +
  geom_line(data = df52, aes(x = 1:nrow(df52), y = Ca_b2, color = "55b")) +
  scale_color_manual(values = c("55a" = "blue", "55b" = "green")) +
  theme_minimal() +
  labs(x = "time", y = "Excess calcium", color = "Attractor state")

plot51_MLC <- ggplot() +
  geom_line(data = df51, aes(x = 1:nrow(df51), y = MLC, color = "55a")) +
  geom_line(data = df52, aes(x = 1:nrow(df52), y = MLC, color = "55b")) +
  scale_color_manual(values = c("55a" = "blue", "55b" = "green")) +
  theme_minimal() +
  labs(x = "time", y = "MLC activity", color = "Attractor state")

patch <- (plot51_AC + plot51_JPs) / (plot51_NO + plot51_Ca1) / (plot51_Ca2 + plot51_MLC)  

# Define where to save the plot for attractor 51

a55_path <- paste('Figures\\ode\\', as.character(55), '_dynamic.pdf', sep='')

# Open the PDF device
pdf(a55_path, width = 11, height = 8)

# Generate the plots
print(patch)

# Close the device to save the file
dev.off()

#####################
# Plot for attractor 62

df75 <- readRDS("Figures\\ode\\75.rds")
df76 <- readRDS("Figures\\ode\\76.rds")
df77 <- readRDS("Figures\\ode\\77.rds")
df78 <- readRDS("Figures\\ode\\78.rds")

plot62_AC <- ggplot() +
  geom_line(data = df75, aes(x = 1:nrow(df75), y = AC, color = "61a")) +
  geom_line(data = df76, aes(x = 1:nrow(df76), y = AC, color = "61b")) +
  geom_line(data = df77, aes(x = 1:nrow(df77), y = AC, color = "61c")) +
  geom_line(data = df78, aes(x = 1:nrow(df78), y = AC, color = "61d")) +
  scale_color_manual(values = c("61a" = "blue", "61b" = "green", "61c" = "red", "61d" = "purple")) +
  theme_minimal() +
  labs(x = "time", y = "AC activity", color = "Attractor state")

plot62_JPs <- ggplot() +
  geom_line(data = df75, aes(x = 1:nrow(df75), y = JPs, color = "61a")) +
  geom_line(data = df76, aes(x = 1:nrow(df76), y = JPs, color = "61b")) +
  geom_line(data = df77, aes(x = 1:nrow(df77), y = JPs, color = "61c")) +
  geom_line(data = df78, aes(x = 1:nrow(df78), y = JPs, color = "61d")) +
  scale_color_manual(values = c("61a" = "blue", "61b" = "green", "61c" = "red", "61d" = "purple")) +
  theme_minimal() +
  labs(x = "time", y = "JPs function", color = "Attractor state")

plot62_RAC1 <- ggplot() +
  geom_line(data = df75, aes(x = 1:nrow(df75), y = RAC1, color = "61a")) +
  geom_line(data = df76, aes(x = 1:nrow(df76), y = RAC1, color = "61b")) +
  geom_line(data = df77, aes(x = 1:nrow(df77), y = RAC1, color = "61c")) +
  geom_line(data = df78, aes(x = 1:nrow(df78), y = RAC1, color = "61d")) +
  scale_color_manual(values = c("61a" = "blue", "61b" = "green", "61c" = "red", "61d" = "purple")) +
  theme_minimal() +
  labs(x = "time", y = "RAC1 activity", color = "Attractor state")

plot62_NO_b1 <- ggplot() +
  geom_line(data = df75, aes(x = 1:nrow(df75), y = NO_b1, color = "61a")) +
  geom_line(data = df76, aes(x = 1:nrow(df76), y = NO_b1, color = "61b")) +
  geom_line(data = df77, aes(x = 1:nrow(df77), y = NO_b1, color = "61c")) +
  geom_line(data = df78, aes(x = 1:nrow(df78), y = NO_b1, color = "61d")) +
  scale_color_manual(values = c("61a" = "blue", "61b" = "green", "61c" = "red", "61d" = "purple")) +
  theme_minimal() +
  labs(x = "time", y = "Normal NO production", color = "Attractor state")

plot62_NO_b2 <- ggplot() +
  geom_line(data = df75, aes(x = 1:nrow(df75), y = NO_b2, color = "61a")) +
  geom_line(data = df76, aes(x = 1:nrow(df76), y = NO_b2, color = "61b")) +
  geom_line(data = df77, aes(x = 1:nrow(df77), y = NO_b2, color = "61c")) +
  geom_line(data = df78, aes(x = 1:nrow(df78), y = NO_b2, color = "61d")) +
  scale_color_manual(values = c("61a" = "blue", "61b" = "green", "61c" = "red", "61d" = "purple")) +
  theme_minimal() +
  labs(x = "time", y = "Excess NO production", color = "Attractor state")

plot62_Ca_b1 <- ggplot() +
  geom_line(data = df75, aes(x = 1:nrow(df75), y = Ca_b1, color = "61a")) +
  geom_line(data = df76, aes(x = 1:nrow(df76), y = Ca_b1, color = "61b")) +
  geom_line(data = df77, aes(x = 1:nrow(df77), y = Ca_b1, color = "61c")) +
  geom_line(data = df78, aes(x = 1:nrow(df78), y = Ca_b1, color = "61d")) +
  scale_color_manual(values = c("61a" = "blue", "61b" = "green", "61c" = "red", "61d" = "purple")) +
  theme_minimal() +
  labs(x = "time", y = "Normal calcium production", color = "Attractor state")

plot62_Ca_b2 <- ggplot() +
  geom_line(data = df75, aes(x = 1:nrow(df75), y = Ca_b2, color = "61a")) +
  geom_line(data = df76, aes(x = 1:nrow(df76), y = Ca_b2, color = "61b")) +
  geom_line(data = df77, aes(x = 1:nrow(df77), y = Ca_b2, color = "61c")) +
  geom_line(data = df78, aes(x = 1:nrow(df78), y = Ca_b2, color = "61d")) +
  scale_color_manual(values = c("61a" = "blue", "61b" = "green", "61c" = "red", "61d" = "purple")) +
  theme_minimal() +
  labs(x = "time", y = "Excess calcium production", color = "Attractor state")

plot62_RHOA <- ggplot() +
  geom_line(data = df75, aes(x = 1:nrow(df75), y = RHOA, color = "61a")) +
  geom_line(data = df76, aes(x = 1:nrow(df76), y = RHOA, color = "61b")) +
  geom_line(data = df77, aes(x = 1:nrow(df77), y = RHOA, color = "61c")) +
  geom_line(data = df78, aes(x = 1:nrow(df78), y = RHOA, color = "61d")) +
  scale_color_manual(values = c("61a" = "blue", "61b" = "green", "61c" = "red", "61d" = "purple")) +
  theme_minimal() +
  labs(x = "time", y = "RHOA activity", color = "Attractor state")

plot62_MLC <- ggplot() +
  geom_line(data = df75, aes(x = 1:nrow(df75), y = MLC, color = "61a")) +
  geom_line(data = df76, aes(x = 1:nrow(df76), y = MLC, color = "61b")) +
  geom_line(data = df77, aes(x = 1:nrow(df77), y = MLC, color = "61c")) +
  geom_line(data = df78, aes(x = 1:nrow(df78), y = MLC, color = "61d")) +
  scale_color_manual(values = c("61a" = "blue", "61b" = "green", "61c" = "red", "61d" = "purple")) +
  theme_minimal() +
  labs(x = "time", y = "MLC activity", color = "Attractor state")

patch1 <- (plot62_AC + plot62_JPs + plot62_RAC1) / (plot62_NO_b1 + plot62_NO_b2 + plot62_Ca_b1) / (plot62_Ca_b2 + plot62_RHOA + plot62_MLC)

a61_path <- paste('Figures\\ode\\', as.character(61), '_dynamic.pdf', sep='')

# Open the PDF device
pdf(a61_path, width = 13, height = 7)

# Generate the plots
print(patch1)

# Close the device to save the file
dev.off()

df55 <- readRDS("Figures\\ode\\55.rds")
df56 <- readRDS("Figures\\ode\\56.rds")

plot59_RHOA <- ggplot() +
  geom_line(data = df55, aes(x = 1:nrow(df55), y = RHOA, color = "59a")) +
  geom_line(data = df56, aes(x = 1:nrow(df56), y = RHOA, color = "59b")) +
  scale_color_manual(values = c("59a" = "blue", "59b" = "green")) +
  theme_minimal() +
  labs(x = "time", y = "RHOA activity", color = "Attractor state")

plot59_RAC1 <- ggplot() +
  geom_line(data = df55, aes(x = 1:nrow(df55), y = RAC1, color = "59a")) +
  geom_line(data = df56, aes(x = 1:nrow(df56), y = RAC1, color = "59b")) +
  scale_color_manual(values = c("59a" = "blue", "59b" = "green")) +
  theme_minimal() +
  labs(x = "time", y = "RAC1 activity", color = "Attractor state")

plot59_MLC <- ggplot() +
  geom_line(data = df55, aes(x = 1:nrow(df55), y = MLC, color = "59a")) +
  geom_line(data = df56, aes(x = 1:nrow(df56), y = MLC, color = "59b")) +
  scale_color_manual(values = c("59a" = "blue", "59b" = "green")) +
  theme_minimal() +
  labs(x = "time", y = "MLC activity", color = "Attractor state")

a59_path <- paste('Figures\\ode\\', as.character(59), '_dynamic.pdf', sep='')
patch2 <- plot59_RHOA + plot59_RAC1 + plot59_MLC
# Open the PDF device
pdf(a59_path, width = 13, height = 5)

# Generate the plots
print(patch2)

# Close the device to save the file
dev.off()