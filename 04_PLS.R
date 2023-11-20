# ===========================================================================================
# Title: Afrotheria - Principal component analysis (PCA)
# Author: Anne Le Maitre, Nicole D. S. Grunstra
# Date: 2023-10-23
# -------------------------------------------------------------------------------------------
# R version: 4.3.1
# Required extensions: ape, geomorph, Morpho, ggrepel, phytools
# -------------------------------------------------------------------------------------------
# Aims: Perform the two-block partial least squares (2B-PLS) analysis
# of the Procrustes-aligned coordinates of 124 landmarks placed on the bony labyrinth 
# against 12 contextual variables for 40 taxa (20 Afrotheria, 20 non-afrotherian), 
# as well as cross-validation of the results, and the phylogenetic PLS analyses
# Input: 
# - Results of the script "00_Definitions.R"
# - Slid landmark coordinates as a CSV file
# - Phylogenetic tree for the 40 taxa
# Output: 
# - 2B-PLS results and related plots
# - Leave-one-out and leave-two-out cross-validation of the PLS loadings
# - Phylogenetic signal in the PLS scores
# - phylogenetic PLS results and related plots
# ===========================================================================================

# ---- Packages ----

library(ape)  # for phylogenetic analyses
library(geomorph)  # geometric morphometrics
library(Morpho)  # geometric morphometrics
library(phytools)  # phylogenetic signal
library(ggrepel)  # plots

# ---- 1.3 Prepare contextual data ----

# Center and scale contextual variables
ECOL <- as.matrix(Afrotheria[, -c(1:4)])  # matrix of the 12 contextual variables
rownames(ECOL) <- Afrotheria$BinomialName
ECOL_std <- scale(ECOL)  # center and scale

# Check associations among contextual variables
PCA_ecol <- prcomp(ECOL_std)
summary(PCA_ecol)
screeplot(PCA_ecol, xlab = "Dimensions", las = 1)  # scree plot
biplot(PCA_ecol, c(1, 2))  # biplot of PC 1 and 2

# ---- 1.4 Prepare landmark data ----

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

# Generalized Procrustes superimposition (no sliding needed)
A.gpa <- procSym(A_slid)

# ---- 1.5 Phylogenetic tree ----

# Load tree
tree <- read.nexus("Grunstra_et_al_Afrotheria_phyl_tree.txt")

# Plot tree
plot(tree, cex = .8)

# ---- 2. Two-block PLS analysis ----

# ---- 2.1 Perform the analysis ----

# 2B-PLS
afro.pls <- pls2B(A.gpa$rotated, ECOL_std)  # analysis
rownames(afro.pls$svd$v) <- colnames(ECOL_std)  # contextual variable names

# Table of the proportion of covariance explained by each PLS dimension
afro.pls$CoVar

# Scree plot
plot(afro.pls$CoVar[, 2], type = "b", col = "blue", las = 1, 
     main = "Scree plot", xlab = "Dimension", 
     ylab = "% of total covariance explained")

# ---- 2.2 Scatter plot: PLS scores ----

# PLS scores for a given dimension
pls <- 1  # example: PLS 1
PLSscores <- cbind(afro.pls$Xscores[, pls], afro.pls$Yscores[, pls])
colnames(PLSscores) <- c("Xscores", "Yscores")
PLSscores <- data.frame(PLSscores)
PLSscores$Habitat <- ecol
PLSscores$Clade <- afro

# Scatter plot of shape scores vs. context scores
ggplot(data = PLSscores, aes(x = Xscores, 
                             y = Yscores, 
                             shape = Clade, colour = Habitat, 
                             label = Afrotheria$CommonName)) +
  geom_text_repel(box.padding = 0.5,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20, 
                  color = "black", 
                  max.overlaps = Inf) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("darkgreen", "orange", "lightblue", "darkblue", "darkred")) +
  theme_classic() + 
  xlab("shape scores") +
  ylab("context scores")

# ---- 2.3 Bar plot: PLS loadings of contextual variables ----

# Contextual variables
plsi <- 1  # dimension
par(mar = c(6, 4, 4, 2) + 0.1)  # set margin sizes
barplot(afro.pls$svd$v[, plsi], 
        main = paste("PLS dimension", plsi), 
        col = "lightblue", las = 2, 
        cex.names = 0.8)
# Back to initial margin values
par(mar = c(5, 4, 4, 2) + 0.1)

# ---- 2.4 3D plots: PLS loadings of Procrustes coordinates ----

# Option 1: use the Morpho function

# Compute coordinates to visualize shape patterns
plsi <- 1  # PLS 1
std <- 2  # 2 standard deviations of the scores
pls.change <- plsCoVar(afro.pls, i = plsi, sdx = std)  

# Visualize 3D shape patterns
layout3d(matrix(1:2, ncol = 2, byrow = TRUE))  # 2 plots in the same window
plotLaby(pls.change$x[,, 1], title = paste("-", std, "sd"))  # Negative scores
next3d()  # Go to the next 3D plot
plotLaby(pls.change$x[,, 2], title = paste("+", std, "sd"))  # Positive scores

# Option 2: use the geomorph function

# Do the 2B-PLS with the geomorph function (same results, except sign changes)
afro.pls2 <- two.b.pls(A.gpa$rotated, ECOL_std)

# Use the function plot.pls to compute required parameters
P <- plot(afro.pls2)
preds <- shape.predictor(P$A1, 
                         x = P$plot.args$x,
                         min = min(P$plot.args$x), 
                         max = max(P$plot.args$x))

# Visualize 3D shape patterns
open3d()
layout3d(matrix(1:2, ncol = 2, byrow = TRUE))  # 2 plots in the same window
plotLaby(preds$min, title = paste("minimum"))  # minimum values
next3d()  # next plot
plotLaby(preds$max, title = paste("maximum"))  # maximum values

# ---- 3. Cross-validation of the PLS scores ----

# ---- 3.1 Leave-one-out cross-validation ----

# Here we perform a cross-validation of the PLS loadings for contextual variables, 
# by leaving out one species at a time

# PLS on the whole sample
afro.pls <- pls2B(A.gpa$rotated, ECOL_std)
rownames(afro.pls$svd$v) <- colnames(ECOL_std)  # contextual variable names

# Ecological variables
plsi <- 1  # dimension
par(mar = c(7, 4, 4, 2) + 0.1)  # set margin sizes
barplot(afro.pls$svd$v[, plsi], 
        main = paste("PLS dimension", plsi), 
        col = "lightblue", las = 2, 
        ylim = c(-0.7, 0.7),
        cex.names = 1)

# PLS with leave-one-out
# (specimen removed *after* GPA and standardization of contextual variables)
for (i in 1:nrow(ECOL_std)) {
  afro.plsi <- pls2B(A.gpa$rotated[, , -i], ECOL_std[-i, ])
  ecoli <- afro.plsi$svd$v
  si <- cor(afro.pls$svd$v[, plsi], afro.plsi$svd$v[, plsi])
  if (si < 0) {ecoli <- -ecoli}
  points(x = 1:nrow(ecoli) * 1.2 - 0.5, y = ecoli[, plsi], 
         pch = 16, col = "blue")
}

# ---- 3.2 Leave-two-out cross-validation ----

# Here we perform a cross-validation of the PLS loadings for contextual variables, 
# by leaving out two species at a time

# PLS on the whole sample
afro.pls <- pls2B(A.gpa$rotated, ECOL_std)
rownames(afro.pls$svd$v) <- colnames(ECOL_std)  # contextual variable names

# PLS with leave-two-out
plsi <- 1  # dimension
par(mar = c(7, 4, 4, 2) + 0.1)  # set margin sizes
barplot(afro.pls$svd$v[, plsi], 
        main = paste("PLS dimension", plsi), 
        col = "lightblue", las = 2, 
        ylim = c(-0.7, 0.7),
        cex.names = 1)

# PLS with leave-two-out
# (specimen removed *after* GPA and standardization of contextual variables)
for (i in 2:nrow(ECOL_std)) {
  for (j in 1:(i-1)) {
    spi <- c(i,j)
    afro.plsi <- pls2B(A.gpa$rotated[, , -spi], ECOL_std[-spi, ])
    ecoli <- afro.plsi$svd$v
    si <- cor(afro.pls$svd$v[, plsi], afro.plsi$svd$v[, plsi])
    if (si < 0) {ecoli <- -ecoli}
    points(x = 1:nrow(ecoli) * 1.2 - 0.5, y = ecoli[, plsi], 
           pch = 16, col = "blue")
  }
}

# Back to initial margin values
par(mar = c(5, 4, 4, 2) + 0.1)

# ---- 4. Phylogenetic signal on PLS scores (phylogeny-unadjusted) ----

# Here we compute the phylogenetic signal in the first four dimensions of the 2B-PLS

# Function to generate phylogenetic signal for multiple traits simultaneously,
# using phytools:phylosig (code copied from Liam Revell's blog)

# Choose the first 4 PLS dimensions
plsi <- 1:4

# Shape scores from 2B-PLS (not shown)
K_shape <- t(simplify2array(apply(afro.pls$Xscores[, plsi], 2, 
                                  phylosig, 
                                  tree = tree, 
                                  method = "K", 
                                  test = TRUE, 
                                  nsim = 1000)))
K_shape
# SW 1: K=0.518 (p=0.043), SW 2: K=0.704 (p=0.004), SW 3: K=0.476 (p=0.023), SW 4: K=0.433 (p=0.052)

# Context scores from same 2B-PLS as above (not shown)
K_context <- t(simplify2array(apply(afro.pls$Yscores[, plsi], 2, 
                                    phylosig, 
                                    tree=tree, 
                                    method="K", 
                                    test=TRUE, 
                                    nsim=1000)))
K_context
# along PLS 1: K=0.689 (p=0.007), PLS 2: K=0.486 (p=0.016), PLS 3: K=0.481 (p=0.017), PLS 4: K=0.385 (p=0.09)

# ---- 5. Phylogenetic PLS ----

# ---- 5.1 Perform the analysis ----

# Phylogenetic PLS
afro.phyl.pls <- phylo.integration(A.gpa$rotated, ECOL_std, phy = tree)
rownames(afro.phyl.pls$svd$v) <- colnames(ECOL_std)  # contextual variable names

# Summary results
summary(afro.phyl.pls)

# Proportion of covariance explained by each PLS dimension
afro.phyl.pls$CoVar <- cbind(afro.phyl.pls$svd$d, 
                             100 * afro.phyl.pls$svd$d^2 / sum(afro.phyl.pls$svd$d^2))
dimnames(afro.phyl.pls$CoVar) <- list(1:12, c("singular value", "% total covar."))

# Scree plot
plot(afro.phyl.pls$CoVar[, 2], type = "b", col = "blue", las = 1, 
     main = "Scree plot", xlab = "Dimension", 
     ylab = "% of total covariance explained")

# ---- 5.2 Scatter plot: phylogenetic PLS scores ----

# PLS scores for a given dimension
pls <- 1  # example: PLS 1
PLSscores <- cbind(afro.phyl.pls$XScores[, pls], afro.phyl.pls$YScores[, pls])
colnames(PLSscores) <- c("Xscores", "Yscores")
PLSscores <- data.frame(PLSscores)
PLSscores$Habitat <- ecol
PLSscores$Clade <- afro

# Scatter plot of shape scores vs. context scores
ggplot(data = PLSscores, aes(x = Xscores, 
                             y = Yscores, 
                             shape = Clade, colour = Habitat, 
                             label = Afrotheria$CommonName)) +
  geom_text_repel(box.padding = 0.5,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20, 
                  color = "black", 
                  max.overlaps = Inf) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("darkgreen", "orange", "lightblue", "darkblue", "darkred")) +
  theme_classic() + 
  xlab("shape scores") +
  ylab("context scores")

# ---- 5.3 Bar plot: phylogenetic PLS loadings of contextual variables ----

# Ecological variables
plsi <- 1  # dimension
par(mar = c(6, 4, 4, 2) + 0.1)  # set margin sizes
barplot(afro.phyl.pls$right.pls.vectors[, plsi], 
        main = paste("PLS dimension", plsi), 
        col = "lightblue", las = 2, 
        cex.names = 0.8)
# Back to initial margin values
par(mar = c(5, 4, 4, 2) + 0.1)

# ---- 5.4 3D plots: phylogenetic PLS loadings of Procrustes coordinates ----

# Compute shape patterns (PLS 1)
phyl.P <- plot(afro.phyl.pls)
phyl.preds <- shape.predictor(phyl.P$A1, 
                              x = phyl.P$plot.args$x,
                              min = min(phyl.P$plot.args$x), 
                              max = max(phyl.P$plot.args$x))

# Visualize 3D shape patterns
layout3d(matrix(1:2, ncol = 2, byrow = TRUE))  # 2 plots in the same window
plotLaby(phyl.preds$min, title = paste("minimum"))  # minimum values
next3d()  # next plot
plotLaby(phyl.preds$max, title = paste("maximum"))  # maximum values

