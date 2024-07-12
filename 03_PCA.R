# ===========================================================================================
# Title: Afrotheria - Principal component analysis (PCA)
# Author: Anne Le Maitre
# Date: 2024-07-12
# -------------------------------------------------------------------------------------------
# R version: 4.3.1
# Required extensions: geomorph, Morpho, ggrepel
# -------------------------------------------------------------------------------------------
# Aims: Perform the principal component analysis (PCA) of the Procrustes-aligned coordinates
# of the 13 anatomical landmarks and 111 semilandmarks placed on the bony labyrinth 
# for 40 taxa (20 Afrotheria, 20 non-afrotherian)
# Input: 
# - Results of the script "00_Definitions.R"
# - Slid landmark coordinates as a CSV file
# Output: PCA results and related plots
# ===========================================================================================

# ---- Packages ----

library(geomorph)  # geometric morphometrics
library(Morpho)  # geometric morphometrics
library(ggrepel)  # plots

# ---- Load landmark data ----

# First, run the R script '00_Definitions'
source("00_Definitions.R")

# Here we directly load the slid landmark coordinates
# They are in the original coordinate space

# Load the slid landmark coordinates
slid_data <- read.csv("Grunstra_et_al_Afrotheria_landmarks_slid_not_aligned.csv")

# Conversion to a matrix
slid <- as.matrix(slid_data[, -1])

# Name rows and columns
rownames(slid) <- slid_data$Species

# Conversion to an array
A_slid <- arrayspecs(slid, p = ie$p, k = 3)

# ---- GPA and PCA ----

# Generalized Procrustes superimposition (no sliding needed) + PCA
A.gpa <- procSym(A_slid)

# Table of the proportion of variance explained by the 10 first PCs
head(A.gpa$Variance, n = 10)

# Scree plot (20 first dimensions)
plot(A.gpa$Variance[1:20, 2], type = "b", col = "blue", las = 1, 
     main = "Scree plot", xlab = "Dimension", 
     ylab = "% of total variance explained")

# ---- Scatter plot of the PC scores ----

# Preparation
variance <- round(A.gpa$Variance[, 2], digits = 1)
PCscores <- data.frame(A.gpa$PCscores)
PCscores$Habitat <- ecol
PCscores$Clade <- afro

# PC1 vs. PC2
ggplot(data = PCscores, aes(x = PC1, 
                            y = PC2, 
                            shape = Clade, colour = Habitat, 
                            label = Afrotheria$CommonName)) +
  coord_fixed(ratio = 1) +
  geom_text_repel(box.padding = 0.5,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20, 
                  color = "black", 
                  max.overlaps = Inf) +
  geom_point(size = 2.5) +
  scale_color_manual(values = couleurs) +
  theme_classic() + 
  xlab(paste("PC1 (", variance[1], "%)", sep = "")) +
  ylab(paste("PC2 (", variance[2], "%)", sep = ""))

# PC3 vs. PC4
ggplot(data = PCscores, aes(x = PC3, 
                            y = PC4, 
                            shape = Clade, colour = Habitat, 
                            label = Afrotheria$CommonName)) +
  coord_fixed(ratio = 1) +
  geom_text_repel(box.padding = 0.5,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20, 
                  color = "black", 
                  max.overlaps = Inf) +
  geom_point(size = 2.5) +
  scale_color_manual(values = couleurs) +
  theme_classic() + 
  xlab(paste("PC3 (", variance[3], "%)", sep = "")) +
  ylab(paste("PC4 (", variance[4], "%)", sep = ""))

# ---- Shape patterns corresponding to each PC ----

# Compute coordinates to visualize shape patterns

# Shape changes (example for - or + 2 standard deviations)
A <- array(NA, dim = c(ie$p, 3, 2))  # prepare an array
pc <- 1  # for PC1
std <- 2  # compute for + or - 2 standard deviations (std)

# Pattern for - 2 std
A[, , 1] <- restoreShapes(scores = - std * sd(A.gpa$PCscores[,pc]), 
                          PC = A.gpa$PCs[, pc], 
                          mshape = A.gpa$mshape)

# Pattern for + 2 std
A[, , 2] <- restoreShapes(scores = std * sd(A.gpa$PCscores[,pc]), 
                          PC = A.gpa$PCs[, pc], 
                          mshape = A.gpa$mshape)

# Visualize 3D shape patterns
layout3d(matrix(1:2, ncol = 2, byrow = TRUE))  # 2 plots in the same window
plotLaby(A[, , 1], title = paste("-", std, "sd"))  # Negative scores
next3d()  # Go to the next 3D plot
plotLaby(A[, , 2], title = paste("+", std, "sd"))  # Positive scores

# ---- OPTIONAL Export data ----

#write.csv(A.gpa$PCscores, file = "Grunstra_et_al_Afrotheria_PCA_Scores.csv")
#write.csv(A.gpa$PCs, file = "Grunstra_et_al_Afrotheria_PCA_Loadings.csv")

